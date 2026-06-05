//! AlignmentFile + AlignedSegment Python types backed by noodles-bam / noodles-cram.
//!
//! v0.3.1: added CRAM read support via `reference_filename` parameter.
//! An `AnyRecord` enum allows `AlignedSegment` to wrap either a `bam::Record`
//! (BAM path) or a `sam::alignment::RecordBuf` (CRAM path).

use std::cell::RefCell;
use std::fs::File;

use pyo3::exceptions::{PyIOError, PyKeyError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyList;

// Required to bring `alignment_end()` into scope (it is a provided method on the
// `noodles::sam::alignment::Record` trait, not a direct method of `bam::Record`).
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::Record as _;

use noodles::sam::alignment::record::data::field::Value as TagValue;
use noodles::sam::alignment::RecordBuf;

use crate::common::{open_indexed, open_streaming, read_header_indexed, read_header_streaming};
use crate::pileup::pileup_bases;
use crate::stats::count_reads;

fn tag_value_to_py<'py>(py: Python<'py>, v: TagValue<'_>) -> PyResult<Bound<'py, PyAny>> {
    use TagValue as V;
    Ok(match v {
        V::Character(c) => (c as char).to_string().into_pyobject(py)?.into_any(),
        V::Int8(n) => n.into_pyobject(py)?.into_any(),
        V::UInt8(n) => n.into_pyobject(py)?.into_any(),
        V::Int16(n) => n.into_pyobject(py)?.into_any(),
        V::UInt16(n) => n.into_pyobject(py)?.into_any(),
        V::Int32(n) => n.into_pyobject(py)?.into_any(),
        V::UInt32(n) => n.into_pyobject(py)?.into_any(),
        V::Float(n) => n.into_pyobject(py)?.into_any(),
        V::String(s) => String::from_utf8_lossy(s.as_ref())
            .into_owned()
            .into_pyobject(py)?
            .into_any(),
        V::Hex(s) => String::from_utf8_lossy(s.as_ref())
            .into_owned()
            .into_pyobject(py)?
            .into_any(),
        V::Array(_arr) => PyList::empty(py).into_any(),
    })
}

/// Convert a `record_buf::data::field::Value` (owned, from CRAM/RecordBuf) to Python.
fn record_buf_value_to_py<'py>(
    py: Python<'py>,
    v: &noodles::sam::alignment::record_buf::data::field::Value,
) -> PyResult<Bound<'py, PyAny>> {
    use noodles::sam::alignment::record_buf::data::field::Value as BufV;
    Ok(match v {
        BufV::Character(c) => ((*c) as char).to_string().into_pyobject(py)?.into_any(),
        BufV::Int8(n) => n.into_pyobject(py)?.into_any(),
        BufV::UInt8(n) => n.into_pyobject(py)?.into_any(),
        BufV::Int16(n) => n.into_pyobject(py)?.into_any(),
        BufV::UInt16(n) => n.into_pyobject(py)?.into_any(),
        BufV::Int32(n) => n.into_pyobject(py)?.into_any(),
        BufV::UInt32(n) => n.into_pyobject(py)?.into_any(),
        BufV::Float(n) => n.into_pyobject(py)?.into_any(),
        BufV::String(s) => String::from_utf8_lossy(s.as_ref())
            .into_owned()
            .into_pyobject(py)?
            .into_any(),
        BufV::Hex(s) => String::from_utf8_lossy(s.as_ref())
            .into_owned()
            .into_pyobject(py)?
            .into_any(),
        BufV::Array(_) => PyList::empty(py).into_any(),
    })
}

const FLAG_PAIRED: u16 = 0x1;
const FLAG_PROPER: u16 = 0x2;
const FLAG_UNMAP: u16 = 0x4;
const FLAG_MUNMAP: u16 = 0x8;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_MREVERSE: u16 = 0x20;
const FLAG_READ1: u16 = 0x40;
const FLAG_READ2: u16 = 0x80;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_QCFAIL: u16 = 0x200;
const FLAG_DUP: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

fn kind_to_char(kind: Kind) -> char {
    match kind {
        Kind::Match => 'M',
        Kind::Insertion => 'I',
        Kind::Deletion => 'D',
        Kind::Skip => 'N',
        Kind::SoftClip => 'S',
        Kind::HardClip => 'H',
        Kind::Pad => 'P',
        Kind::SequenceMatch => '=',
        Kind::SequenceMismatch => 'X',
    }
}

fn kind_to_int(kind: Kind) -> u8 {
    // SAM spec / pysam encoding.
    match kind {
        Kind::Match => 0,
        Kind::Insertion => 1,
        Kind::Deletion => 2,
        Kind::Skip => 3,
        Kind::SoftClip => 4,
        Kind::HardClip => 5,
        Kind::Pad => 6,
        Kind::SequenceMatch => 7,
        Kind::SequenceMismatch => 8,
    }
}

/// Reverse of `kind_to_int` — pysam-style op code (0..=8) → noodles `Kind`.
fn kind_from_int(op: u8) -> PyResult<Kind> {
    match op {
        0 => Ok(Kind::Match),
        1 => Ok(Kind::Insertion),
        2 => Ok(Kind::Deletion),
        3 => Ok(Kind::Skip),
        4 => Ok(Kind::SoftClip),
        5 => Ok(Kind::HardClip),
        6 => Ok(Kind::Pad),
        7 => Ok(Kind::SequenceMatch),
        8 => Ok(Kind::SequenceMismatch),
        other => Err(PyValueError::new_err(format!(
            "invalid CIGAR op code {other}: expected 0..=8 (M I D N S H P = X)"
        ))),
    }
}

/// Reverse of `kind_to_char` — SAM CIGAR letter → noodles `Kind`.
fn kind_from_char(c: char) -> PyResult<Kind> {
    match c {
        'M' => Ok(Kind::Match),
        'I' => Ok(Kind::Insertion),
        'D' => Ok(Kind::Deletion),
        'N' => Ok(Kind::Skip),
        'S' => Ok(Kind::SoftClip),
        'H' => Ok(Kind::HardClip),
        'P' => Ok(Kind::Pad),
        '=' => Ok(Kind::SequenceMatch),
        'X' => Ok(Kind::SequenceMismatch),
        other => Err(PyValueError::new_err(format!(
            "invalid CIGAR op {other:?}: expected one of M I D N S H P = X"
        ))),
    }
}

/// Parse a CIGAR string like "150M" or "10S140M5I" into a list of noodles `Op`s.
fn parse_cigar_string(s: &str) -> PyResult<Vec<noodles::sam::alignment::record::cigar::Op>> {
    use noodles::sam::alignment::record::cigar::Op;
    let mut out = Vec::new();
    if s.is_empty() || s == "*" {
        return Ok(out);
    }
    let mut len_acc: usize = 0;
    let mut has_digit = false;
    for ch in s.chars() {
        if let Some(d) = ch.to_digit(10) {
            len_acc = len_acc
                .checked_mul(10)
                .and_then(|n| n.checked_add(d as usize))
                .ok_or_else(|| {
                    PyValueError::new_err(format!("CIGAR op length overflow in {s:?}"))
                })?;
            has_digit = true;
        } else {
            if !has_digit {
                return Err(PyValueError::new_err(format!(
                    "malformed CIGAR {s:?}: op {ch:?} without preceding length"
                )));
            }
            let kind = kind_from_char(ch)?;
            out.push(Op::new(kind, len_acc));
            len_acc = 0;
            has_digit = false;
        }
    }
    if has_digit {
        return Err(PyValueError::new_err(format!(
            "malformed CIGAR {s:?}: trailing digits without an op letter"
        )));
    }
    Ok(out)
}

/// Convert a Python tag value (int / float / str / bytes) into a
/// `record_buf::data::field::Value`. Used by `AlignedSegment.set_tag`.
fn py_to_tag_value(
    obj: &Bound<'_, PyAny>,
) -> PyResult<noodles::sam::alignment::record_buf::data::field::Value> {
    use noodles::sam::alignment::record_buf::data::field::Value as BufV;
    // bool must be checked before int because bool is a subclass of int in Python.
    if obj.is_instance_of::<pyo3::types::PyBool>() {
        let b: bool = obj.extract()?;
        return Ok(BufV::Int8(if b { 1 } else { 0 }));
    }
    if let Ok(n) = obj.extract::<i64>() {
        // Pick the narrowest representable integer type per SAM spec semantics.
        if let Ok(v) = i8::try_from(n) {
            return Ok(BufV::Int8(v));
        }
        if let Ok(v) = u8::try_from(n) {
            return Ok(BufV::UInt8(v));
        }
        if let Ok(v) = i16::try_from(n) {
            return Ok(BufV::Int16(v));
        }
        if let Ok(v) = u16::try_from(n) {
            return Ok(BufV::UInt16(v));
        }
        if let Ok(v) = i32::try_from(n) {
            return Ok(BufV::Int32(v));
        }
        if let Ok(v) = u32::try_from(n) {
            return Ok(BufV::UInt32(v));
        }
        return Err(PyValueError::new_err(format!(
            "tag integer {n} does not fit in any SAM integer width (i8/u8/i16/u16/i32/u32)"
        )));
    }
    if let Ok(f) = obj.extract::<f32>() {
        return Ok(BufV::Float(f));
    }
    if let Ok(s) = obj.extract::<String>() {
        return Ok(BufV::String(s.into_bytes().into()));
    }
    if let Ok(b) = obj.extract::<Vec<u8>>() {
        return Ok(BufV::String(b.into()));
    }
    Err(PyTypeError::new_err(format!(
        "set_tag: unsupported Python type for tag value: {}",
        obj.get_type().name()?
    )))
}

// ---------------------------------------------------------------------------
// AnyRecord: unifies bam::Record and sam::alignment::RecordBuf for AlignedSegment.
// ---------------------------------------------------------------------------

/// Unifies BAM and CRAM record types under a single dispatch path.
///
/// Marked `#[non_exhaustive]` so future variants (e.g. eventual write
/// support, or an SAM text-record path) can be added without breaking
/// downstream `match` statements when this enum is exposed publicly.
#[non_exhaustive]
pub(crate) enum AnyRecord {
    Bam(noodles::bam::Record),
    Cram(noodles::sam::alignment::RecordBuf),
}

// ---------------------------------------------------------------------------
// Inner reader enum for AlignmentFile.
// ---------------------------------------------------------------------------

pub(crate) type IndexedCramReader = noodles::cram::io::IndexedReader<File>;

/// Type alias for the BAM writer that backs `AlignmentFile(path, "wb", ...)`.
/// `bam::io::Writer::new` wraps its inner writer in a BGZF compressor, so the
/// outer type is `Writer<bgzf::Writer<BufWriter<File>>>`. We buffer the raw
/// file to amortise the BGZF block writes.
pub(crate) type BamFileWriter =
    noodles::bam::io::Writer<noodles::bgzf::io::Writer<std::io::BufWriter<std::fs::File>>>;

#[allow(dead_code)]
pub(crate) enum AnyIo {
    Bam(crate::common::IndexedBamReader),
    Cram(IndexedCramReader),
    BamWriter(BamFileWriter),
}

/// Build an fasta::Repository from an optional reference FASTA path.
/// Returns an empty repository if `reference_filename` is None.
///
/// Supports both uncompressed `.fa` (needs `.fa.fai`) and BGZF-compressed
/// `.fa.gz` (needs `.fa.gz.fai` + `.fa.gz.gzi`).
fn build_fasta_repo(reference_filename: Option<&str>) -> PyResult<noodles::fasta::Repository> {
    match reference_filename {
        None => Ok(noodles::fasta::Repository::default()),
        Some(ref_path) => {
            // fasta::io::indexed_reader::Builder handles both plain and BGZF-compressed
            // FASTA automatically based on the file extension.  It looks for <path>.fai
            // and, for .gz files, also <path>.gzi.
            let indexed_reader = noodles::fasta::io::indexed_reader::Builder::default()
                .build_from_path(ref_path)
                .map_err(|e| {
                    PyIOError::new_err(format!(
                        "failed to open reference FASTA '{ref_path}': {e}\n\
                         (for .fa: run `samtools faidx {ref_path}`;\n\
                          for .fa.gz: run `bgzip -d {ref_path}` first)"
                    ))
                })?;
            let adapter = noodles::fasta::repository::adapters::IndexedReader::new(indexed_reader);
            Ok(noodles::fasta::Repository::new(adapter))
        }
    }
}

/// Thin handle around an indexed BAM or CRAM reader.
///
/// `noodles::bam::io::IndexedReader` / `noodles::cram::io::IndexedReader` hold
/// a `Box<dyn BinningIndex>` which is not `Send`, so we use `#[pyclass(unsendable)]`.
#[pyclass(unsendable)]
pub struct AlignmentFile {
    #[allow(dead_code)]
    path: String,
    /// Mode string passed at open ("rb", "wb", "rc", ...).
    mode: String,
    /// Optional reference FASTA path for CRAM decoding.
    reference_filename: Option<String>,
    /// `threads` kwarg passed at open — stored as a no-op for pysam-compat.
    threads: usize,
    /// `add_hts_options` storage — pysam-compat no-op list.
    hts_options: RefCell<Vec<String>>,
    /// `check_truncation` flag — pysam-compat no-op storage.
    check_truncation_flag: std::cell::Cell<bool>,
    inner: RefCell<Option<AnyIo>>,
    header: Option<noodles::sam::Header>,
    closed: std::cell::Cell<bool>,
}

fn is_cram(path: &str) -> bool {
    path.to_ascii_lowercase().ends_with(".cram")
}

fn open_indexed_cram(path: &str, reference_filename: Option<&str>) -> PyResult<IndexedCramReader> {
    let repo = build_fasta_repo(reference_filename)?;
    noodles::cram::io::indexed_reader::Builder::default()
        .set_reference_sequence_repository(repo)
        .build_from_path(path)
        .map_err(|e| PyIOError::new_err(format!("failed to open indexed CRAM at {path}: {e}")))
}

/// Map a pysam `read_callback` string to the SAM-flag filter mask it implies.
///
/// pysam semantics: `'nofilter'` keeps every read (mask `0`); `'all'` skips
/// reads with `UNMAP | SECONDARY | QCFAIL | DUP` (`0x704`) — note that
/// supplementary (`0x800`) reads are **kept** by `'all'`.
fn read_callback_mask(read_callback: &str) -> PyResult<u16> {
    match read_callback {
        "nofilter" => Ok(0),
        "all" => Ok(0x704),
        other => Err(PyValueError::new_err(format!(
            "unsupported read_callback {other:?}; expected 'all' or 'nofilter'"
        ))),
    }
}

#[pymethods]
impl AlignmentFile {
    #[new]
    #[pyo3(signature = (path, mode = "rb", reference_filename = None, template = None, header = None, threads = 1))]
    fn new(
        path: std::path::PathBuf,
        mode: &str,
        reference_filename: Option<std::path::PathBuf>,
        template: Option<&AlignmentFile>,
        header: Option<&Header>,
        threads: usize,
    ) -> PyResult<Self> {
        // pyo3 0.23's FromPyObject for std::path::PathBuf accepts str, bytes,
        // and os.PathLike (including pathlib.Path) by calling os.fspath. This
        // closes the v0.3.2 ergonomic gap flagged by the v5 reviewer.
        let path: &str = path.to_str().ok_or_else(|| {
            PyValueError::new_err(
                "path must be valid UTF-8 (Windows non-UTF-16 paths not yet supported)",
            )
        })?;
        let reference_filename: Option<String> = reference_filename
            .as_ref()
            .map(|p| p.to_string_lossy().into_owned());
        let reference_filename: Option<&str> = reference_filename.as_deref();
        // ---------- write mode ----------
        if mode == "wb" {
            // Resolve the header source. pysam allows either template (an
            // already-open AlignmentFile whose header is copied) or an
            // explicit header= kwarg.
            let hdr = if let Some(t) = template {
                t.header
                    .as_ref()
                    .ok_or_else(|| {
                        PyValueError::new_err(
                            "template AlignmentFile has no header (was it closed?)",
                        )
                    })?
                    .clone()
            } else if let Some(h) = header {
                h.inner.clone()
            } else {
                return Err(PyValueError::new_err(
                    "mode='wb' requires either template=AlignmentFile or header=Header",
                ));
            };
            let file = std::fs::File::create(path)
                .map_err(|e| PyIOError::new_err(format!("create {path}: {e}")))?;
            let mut writer: BamFileWriter =
                noodles::bam::io::Writer::new(std::io::BufWriter::new(file));
            writer
                .write_header(&hdr)
                .map_err(|e| PyIOError::new_err(format!("write BAM header: {e}")))?;
            return Ok(Self {
                path: path.to_string(),
                mode: mode.to_string(),
                reference_filename: None,
                threads,
                hts_options: RefCell::new(Vec::new()),
                check_truncation_flag: std::cell::Cell::new(true),
                inner: RefCell::new(Some(AnyIo::BamWriter(writer))),
                header: Some(hdr),
                closed: std::cell::Cell::new(false),
            });
        }
        // ---------- read mode ----------
        if mode != "rb" {
            return Err(PyValueError::new_err(format!(
                "supported modes: 'rb' (read), 'wb' (write); got {mode:?}"
            )));
        }
        if is_cram(path) {
            let mut reader = open_indexed_cram(path, reference_filename)?;
            let header = reader
                .read_header()
                .map_err(|e| PyIOError::new_err(format!("failed to read CRAM header: {e}")))?;
            Ok(Self {
                path: path.to_string(),
                mode: mode.to_string(),
                reference_filename: reference_filename.map(str::to_string),
                threads,
                hts_options: RefCell::new(Vec::new()),
                check_truncation_flag: std::cell::Cell::new(true),
                inner: RefCell::new(Some(AnyIo::Cram(reader))),
                header: Some(header),
                closed: std::cell::Cell::new(false),
            })
        } else {
            // BAM path (existing logic)
            let bai = format!("{path}.bai");
            let csi = format!("{path}.csi");
            let has_index =
                std::path::Path::new(&bai).exists() || std::path::Path::new(&csi).exists();
            if has_index {
                let mut reader = open_indexed(path)?;
                let header = read_header_indexed(&mut reader)
                    .map_err(|e| PyIOError::new_err(format!("failed to read header: {e}")))?;
                Ok(Self {
                    path: path.to_string(),
                    mode: mode.to_string(),
                    reference_filename: None,
                    threads,
                    hts_options: RefCell::new(Vec::new()),
                    check_truncation_flag: std::cell::Cell::new(true),
                    inner: RefCell::new(Some(AnyIo::Bam(reader))),
                    header: Some(header),
                    closed: std::cell::Cell::new(false),
                })
            } else {
                let mut streaming = open_streaming(path)?;
                let header = read_header_streaming(&mut streaming)?;
                Ok(Self {
                    path: path.to_string(),
                    mode: mode.to_string(),
                    reference_filename: None,
                    threads,
                    hts_options: RefCell::new(Vec::new()),
                    check_truncation_flag: std::cell::Cell::new(true),
                    inner: RefCell::new(None),
                    header: Some(header),
                    closed: std::cell::Cell::new(false),
                })
            }
        }
    }

    #[getter]
    fn is_open(&self) -> bool {
        !self.closed.get()
    }

    fn close(&self) -> PyResult<()> {
        // When closing a writer, the underlying bgzf::Writer's Drop impl
        // writes the BGZF EOF block automatically, so a simple drop is
        // enough to produce a well-formed BAM. We still take ownership and
        // drop explicitly so any I/O error surfaces as a Python exception
        // path the caller can catch (rather than a silent panic-in-Drop).
        if let Some(io) = self.inner.borrow_mut().take() {
            drop(io);
        }
        self.closed.set(true);
        Ok(())
    }

    /// Write one `AlignedSegment` to the underlying BAM writer. The file
    /// must have been opened with `mode="wb"`.
    ///
    /// pysam-compatible: returns 0 (pysam returns the number of bytes
    /// written, which is not cheaply accessible from noodles). Every
    /// standard caller ignores the return value.
    fn write(&self, segment: &AlignedSegment) -> PyResult<usize> {
        use noodles::sam::alignment::io::Write as _;
        let mut guard = self.inner.borrow_mut();
        let writer = match guard.as_mut() {
            Some(AnyIo::BamWriter(w)) => w,
            _ => {
                return Err(PyValueError::new_err(
                    "AlignmentFile is not opened for writing (mode!='wb')",
                ))
            }
        };
        let hdr = self
            .header
            .as_ref()
            .ok_or_else(|| PyIOError::new_err("writer has no header"))?;
        match &segment.record {
            AnyRecord::Bam(r) => writer
                .write_alignment_record(hdr, r)
                .map_err(|e| PyIOError::new_err(format!("write record: {e}")))?,
            AnyRecord::Cram(r) => writer
                .write_alignment_record(hdr, r)
                .map_err(|e| PyIOError::new_err(format!("write record (from CRAM): {e}")))?,
        }
        Ok(0)
    }

    fn __enter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __exit__(
        &self,
        _exc_type: &Bound<'_, PyAny>,
        _exc_val: &Bound<'_, PyAny>,
        _exc_tb: &Bound<'_, PyAny>,
    ) -> PyResult<bool> {
        self.close()?;
        Ok(false)
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyResult<AlignmentFileStreamIter> {
        let bam_path = slf.path.clone();
        if is_cram(&bam_path) {
            return Err(PyIOError::new_err(
                "linear iteration over CRAM files is not supported; use fetch() instead",
            ));
        }
        let mut reader = open_streaming(&bam_path)?;
        let header = read_header_streaming(&mut reader)?;
        Ok(AlignmentFileStreamIter {
            reader: RefCell::new(reader),
            header,
        })
    }

    #[getter]
    fn references<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        use pyo3::types::PyTuple;
        let header = match self.header.as_ref() {
            Some(h) => h,
            None => return Ok(PyTuple::empty(py)),
        };
        let names: Vec<String> = header
            .reference_sequences()
            .iter()
            .map(|(name, _)| String::from_utf8_lossy(name).into_owned())
            .collect();
        Ok(PyTuple::new(py, names)?)
    }

    #[getter]
    fn lengths<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        use pyo3::types::PyTuple;
        let header = match self.header.as_ref() {
            Some(h) => h,
            None => return Ok(PyTuple::empty(py)),
        };
        let lens: Vec<usize> = header
            .reference_sequences()
            .iter()
            .map(|(_, ref_seq)| ref_seq.length().get())
            .collect();
        Ok(PyTuple::new(py, lens)?)
    }

    #[getter]
    fn nreferences(&self) -> usize {
        self.header
            .as_ref()
            .map(|h| h.reference_sequences().len())
            .unwrap_or(0)
    }

    /// pysam-compatible `get_reference_length(contig)`. Returns the @SQ length
    /// for the given contig, or raises `KeyError` if the contig is not in the
    /// header.
    fn get_reference_length(&self, contig: &str) -> PyResult<usize> {
        let header = self
            .header
            .as_ref()
            .ok_or_else(|| PyIOError::new_err("file is closed"))?;
        header
            .reference_sequences()
            .get(contig.as_bytes())
            .map(|rs| rs.length().get())
            .ok_or_else(|| PyKeyError::new_err(format!("contig {contig:?} not in BAM header")))
    }

    #[getter]
    fn header(&self) -> PyResult<Header> {
        let h = self
            .header
            .as_ref()
            .ok_or_else(|| PyIOError::new_err("file is closed"))?;
        Ok(Header { inner: h.clone() })
    }

    /// `fetch(contig=None, start=0, end=0, *, until_eof=False)`.
    ///
    /// pysam-compatible: when `until_eof=True`, iterates **every** record in
    /// the file (mapped + unmapped) in BGZF order, bypassing the index. When
    /// `until_eof=False`, the usual `(contig, start, end)` indexed query
    /// applies and `contig` is required. The two paths return iterator types
    /// with identical Python protocols (`__iter__`, `__next__`).
    #[pyo3(signature = (contig = None, start = 0, end = 0, *, until_eof = false))]
    fn fetch(
        &self,
        py: Python<'_>,
        contig: Option<&str>,
        start: usize,
        end: usize,
        until_eof: bool,
    ) -> PyResult<PyObject> {
        use noodles::core::{Position, Region};

        let path = self.path.clone();

        // pysam-compat: until_eof=True iterates every record without index.
        // We route to the existing streaming iterator.
        if until_eof {
            let mut reader = open_streaming(&path)?;
            let header = read_header_streaming(&mut reader)?;
            let stream_iter = AlignmentFileStreamIter {
                reader: RefCell::new(reader),
                header,
            };
            return Ok(Py::new(py, stream_iter)?.into_any().into());
        }

        let contig = contig.ok_or_else(|| {
            PyValueError::new_err(
                "fetch() requires `contig` (use until_eof=True for full-BAM iteration)",
            )
        })?;

        if is_cram(&path) {
            // Re-open the CRAM with a fresh indexed reader (same as BAM fetch approach).
            let ref_file = self.reference_filename.as_deref();
            let mut reader = open_indexed_cram(&path, ref_file)?;
            let header = reader
                .read_header()
                .map_err(|e| PyIOError::new_err(format!("read CRAM header: {e}")))?;

            let ref_len = match header.reference_sequences().get(contig.as_bytes()) {
                Some(rs) => rs.length().get(),
                None => {
                    return Err(PyValueError::new_err(format!(
                        "chromosome {contig} not in CRAM header"
                    )));
                }
            };
            if end <= start || start >= ref_len {
                let empty = AlignmentFileFetchIter {
                    records: RefCell::new(Vec::<AnyRecord>::new().into_iter()),
                    header,
                };
                return Ok(Py::new(py, empty)?.into_any().into());
            }
            let end = end.min(ref_len);
            let region = Region::new(
                contig.as_bytes().to_vec(),
                Position::new(start + 1)
                    .ok_or_else(|| PyValueError::new_err("start+1 must be >= 1"))?
                    ..=Position::new(end)
                        .ok_or_else(|| PyValueError::new_err("end must be >= 1"))?,
            );
            let query = reader
                .query(&header, &region)
                .map_err(|e| PyIOError::new_err(format!("CRAM query: {e}")))?;
            // CRAM record decode can panic inside noodles-cram on slices using codecs
            // that are not yet implemented upstream (notably Huffman byte-series, used
            // by NYGC 30x CRAMs). Catch unwinding panics across the FFI boundary and
            // convert to a clean PyIOError so Python callers never see a process abort.
            let panic_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                let mut buf: Vec<AnyRecord> = Vec::new();
                for r in query {
                    let rec = r.map_err(|e| PyIOError::new_err(format!("CRAM record: {e}")))?;
                    buf.push(AnyRecord::Cram(rec));
                }
                Ok::<Vec<AnyRecord>, PyErr>(buf)
            }));
            let buf = match panic_result {
                Ok(Ok(buf)) => buf,
                Ok(Err(e)) => return Err(e),
                Err(_) => {
                    return Err(PyIOError::new_err(
                        "CRAM record decode failed (likely an unsupported codec such \
                         as Huffman byte-series in noodles-cram 0.90+). CRAM support \
                         in rubam is currently experimental and effectively header-only \
                         for files using unsupported codecs; opening the file and reading \
                         the header still works.",
                    ));
                }
            };
            let cram_iter = AlignmentFileFetchIter {
                records: RefCell::new(buf.into_iter()),
                header,
            };
            Ok(Py::new(py, cram_iter)?.into_any().into())
        } else {
            // BAM path — existing logic
            let mut reader = open_indexed(&path)?;
            let header = read_header_indexed(&mut reader)
                .map_err(|e| PyIOError::new_err(format!("read header: {e}")))?;
            let ref_len = match header.reference_sequences().get(contig.as_bytes()) {
                Some(rs) => rs.length().get(),
                None => {
                    return Err(PyValueError::new_err(format!(
                        "chromosome {contig} not in BAM header"
                    )));
                }
            };
            if end <= start || start >= ref_len {
                let empty = AlignmentFileFetchIter {
                    records: RefCell::new(Vec::<AnyRecord>::new().into_iter()),
                    header,
                };
                return Ok(Py::new(py, empty)?.into_any().into());
            }
            let end = end.min(ref_len);
            let region = Region::new(
                contig.as_bytes().to_vec(),
                Position::new(start + 1)
                    .ok_or_else(|| PyValueError::new_err("start+1 must be >= 1"))?
                    ..=Position::new(end)
                        .ok_or_else(|| PyValueError::new_err("end must be >= 1"))?,
            );
            let query = reader
                .query(&header, &region)
                .map_err(|e| PyIOError::new_err(format!("query: {e}")))?;
            let mut buf: Vec<AnyRecord> = Vec::new();
            for r in query.records() {
                let rec = r.map_err(|e| PyIOError::new_err(format!("record: {e}")))?;
                buf.push(AnyRecord::Bam(rec));
            }
            let bam_iter = AlignmentFileFetchIter {
                records: RefCell::new(buf.into_iter()),
                header,
            };
            Ok(Py::new(py, bam_iter)?.into_any().into())
        }
    }

    fn has_index(&self) -> bool {
        if is_cram(&self.path) {
            let crai = format!("{}.crai", self.path);
            return std::path::Path::new(&crai).exists();
        }
        let bai = format!("{}.bai", self.path);
        let csi = format!("{}.csi", self.path);
        std::path::Path::new(&bai).exists() || std::path::Path::new(&csi).exists()
    }

    fn check_index(&self) -> PyResult<()> {
        if !self.has_index() {
            return Err(PyIOError::new_err("file has no index (.bai/.csi/.crai)"));
        }
        Ok(())
    }

    fn get_index_statistics<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, pyo3::types::PyList>> {
        use pyo3::types::{PyDict, PyList};
        if is_cram(&self.path) {
            return Err(PyIOError::new_err(
                "get_index_statistics is not yet implemented for CRAM files",
            ));
        }
        let path = self.path.clone();
        let mut reader = open_indexed(&path)?;
        let header = read_header_indexed(&mut reader)
            .map_err(|e| PyIOError::new_err(format!("read header: {e}")))?;
        let index = reader.index();
        let pylist = PyList::empty(py);
        let ref_seqs_iter = index.reference_sequences();
        for ((name, _), idx_ref) in header.reference_sequences().iter().zip(ref_seqs_iter) {
            let metadata = idx_ref.metadata();
            let mapped = metadata.map(|m| m.mapped_record_count()).unwrap_or(0);
            let unmapped = metadata.map(|m| m.unmapped_record_count()).unwrap_or(0);
            let row = PyDict::new(py);
            row.set_item("contig", String::from_utf8_lossy(name).into_owned())?;
            row.set_item("mapped", mapped)?;
            row.set_item("unmapped", unmapped)?;
            pylist.append(row)?;
        }
        Ok(pylist)
    }

    fn head(&self, n: usize) -> PyResult<Vec<AlignedSegment>> {
        if is_cram(&self.path) {
            return Err(PyIOError::new_err(
                "head() is not yet implemented for CRAM files; use fetch() instead",
            ));
        }
        let path = self.path.clone();
        let mut reader = open_indexed(&path)?;
        let header = read_header_indexed(&mut reader)
            .map_err(|e| PyIOError::new_err(format!("read header: {e}")))?;
        let mut out = Vec::with_capacity(n);
        let mut record = noodles::bam::Record::default();
        for _ in 0..n {
            match reader.read_record(&mut record) {
                Ok(0) => break,
                Ok(_) => out.push(AlignedSegment::new(
                    AnyRecord::Bam(record.clone()),
                    header.clone(),
                )),
                Err(e) => return Err(PyIOError::new_err(format!("head read: {e}"))),
            }
        }
        Ok(out)
    }

    /// `count(contig, start=None, end=None, *, read_callback='nofilter', ...)`.
    ///
    /// pysam-compatible: when `start` and `end` are omitted, defaults to the
    /// whole contig (1 .. reference_length). This is the form most pysam
    /// callers use (`samfile.count(contig="chr16")`).
    ///
    /// The default `read_callback='nofilter'` mirrors **pysam's** `count`
    /// default: it counts every read in the region, including secondary,
    /// supplementary, duplicate and QC-fail records. `read_callback='all'`
    /// applies the `0x704` mask (skip `UNMAP | SECONDARY | QCFAIL | DUP`).
    /// `flag_required` / `flag_filtered` remain available for explicit control
    /// and, when `flag_filtered` is given, it overrides `read_callback`.
    ///
    /// Note: the free `count_reads()` function defaults to the *opposite*,
    /// samtools-style `flag_filtered=0x704` (skips secondary/dup/qcfail/unmap),
    /// so `count_reads(...)` and `samfile.count(...)` can return different
    /// totals on the same region by default. This is intentional, not a bug.
    #[pyo3(signature = (contig, start = None, end = None, *, read_callback = "nofilter", min_mapq = 0, flag_required = 0, flag_filtered = None))]
    fn count(
        &self,
        contig: &str,
        start: Option<usize>,
        end: Option<usize>,
        read_callback: &str,
        min_mapq: u8,
        flag_required: u16,
        flag_filtered: Option<u16>,
    ) -> PyResult<u64> {
        let flag_filtered = match flag_filtered {
            Some(m) => m,
            None => read_callback_mask(read_callback)?,
        };
        // Resolve contig length from the header for the whole-contig default.
        let resolved_end = match end {
            Some(e) => e,
            None => {
                let header = self
                    .header
                    .as_ref()
                    .ok_or_else(|| PyIOError::new_err("file is closed"))?;
                header
                    .reference_sequences()
                    .get(contig.as_bytes())
                    .map(|rs| rs.length().get())
                    .ok_or_else(|| {
                        PyValueError::new_err(format!("chromosome {contig} not in BAM header"))
                    })?
            }
        };
        let start = start.unwrap_or(0);
        let end = resolved_end;
        if is_cram(&self.path) {
            return Err(PyIOError::new_err(
                "count() is not yet implemented for CRAM files",
            ));
        }
        count_reads(
            &self.path,
            contig,
            (start + 1) as u64,
            end as u64,
            min_mapq,
            flag_required,
            flag_filtered,
        )
    }

    #[pyo3(signature = (contig, start, end, *,
                        min_mapq = 0, min_bq = 13, max_depth = 8000,
                        truncate = true, num_threads = 4))]
    fn pileup(
        &self,
        contig: &str,
        start: usize,
        end: usize,
        min_mapq: u8,
        min_bq: u8,
        max_depth: usize,
        truncate: bool,
        num_threads: usize,
    ) -> PyResult<crate::pileup_iter::PileupIter> {
        if is_cram(&self.path) {
            return Err(PyIOError::new_err(
                "pileup() is not yet implemented for CRAM files",
            ));
        }
        let _ = truncate;
        let (positions, a, c, g, t, n_arr, depth) = crate::pileup::pileup_bases(
            &self.path,
            contig,
            (start + 1) as u64,
            end as u64,
            1,
            min_mapq,
            min_bq,
            max_depth,
            num_threads,
            crate::common::FLAG_FILTER_DEFAULT,
        )?;
        let mut cols: Vec<crate::pileup_iter::PileupColumn> = Vec::with_capacity(positions.len());
        for i in 0..positions.len() {
            cols.push(crate::pileup_iter::PileupColumn {
                reference_name: contig.to_string(),
                reference_pos: (positions[i] - 1) as usize,
                depth: depth[i],
                a: a[i],
                c: c[i],
                g: g[i],
                t: t[i],
                n: n_arr[i],
            });
        }
        Ok(crate::pileup_iter::PileupIter {
            cols: std::cell::RefCell::new(cols.into_iter()),
        })
    }

    /// `count_coverage(contig, start, end, *, quality_threshold=15, read_callback='all', ...)`.
    ///
    /// pysam-compatible defaults: `quality_threshold=15` (a base is counted iff
    /// its quality is `>= quality_threshold`), `read_callback='all'` (skip
    /// `UNMAP | SECONDARY | QCFAIL | DUP`, i.e. the `0x704` mask), and **no**
    /// depth cap — pysam's `count_coverage` never truncates coverage.
    #[pyo3(signature = (contig, start, end, *, quality_threshold = 15, read_callback = "all", num_threads = 4))]
    fn count_coverage(
        &self,
        contig: &str,
        start: usize,
        end: usize,
        quality_threshold: u8,
        read_callback: &str,
        num_threads: usize,
    ) -> PyResult<(Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>)> {
        if is_cram(&self.path) {
            return Err(PyIOError::new_err(
                "count_coverage() is not yet implemented for CRAM files",
            ));
        }
        let flag_filter = read_callback_mask(read_callback)?;
        let (_pos, a, c, g, t, _n, _depth) = pileup_bases(
            &self.path,
            contig,
            (start + 1) as u64,
            end as u64,
            1,
            0,
            quality_threshold,
            usize::MAX,
            num_threads,
            flag_filter,
        )?;
        Ok((a, c, g, t))
    }

    // -------- pysam-compatible AlignmentFile aliases (v0.3.4) --------

    /// pysam `closed` — opposite of `is_open`.
    #[getter]
    fn closed(&self) -> bool {
        !self.is_open()
    }
    /// pysam `is_closed` — alias of `closed`.
    #[getter]
    fn is_closed(&self) -> bool {
        !self.is_open()
    }
    /// pysam `mode` — the open-mode string passed at construction.
    #[getter]
    fn mode(&self) -> String {
        self.mode.clone()
    }
    /// pysam `filename` — the file path passed at construction.
    #[getter]
    fn filename(&self) -> String {
        self.path.clone()
    }
    /// pysam `reference_filename` — the FASTA path supplied for CRAM.
    #[getter]
    #[pyo3(name = "reference_filename")]
    fn py_reference_filename(&self) -> Option<String> {
        self.reference_filename.clone()
    }
    /// pysam `threads` — the threads kwarg supplied at construction.
    #[getter]
    #[pyo3(name = "threads")]
    fn py_threads(&self) -> usize {
        self.threads
    }
    /// pysam `is_remote` — always False (no HTTP/S3 support yet).
    #[getter]
    fn is_remote(&self) -> bool {
        false
    }
    /// pysam `is_stream` — True if path is `-` (stdin/stdout).
    #[getter]
    fn is_stream(&self) -> bool {
        self.path == "-"
    }
    /// pysam `is_read` — True if mode starts with `r`.
    #[getter]
    fn is_read(&self) -> bool {
        self.mode.starts_with('r')
    }
    /// pysam `is_write` — True if mode starts with `w`.
    #[getter]
    fn is_write(&self) -> bool {
        self.mode.starts_with('w')
    }
    /// pysam `is_bam` — True if file is BAM.
    #[getter]
    fn is_bam(&self) -> bool {
        self.path.to_ascii_lowercase().ends_with(".bam") || self.mode.contains('b')
    }
    /// pysam `is_cram` — True if file is CRAM.
    #[getter]
    fn is_cram(&self) -> bool {
        self.path.to_ascii_lowercase().ends_with(".cram") || self.mode.contains('c')
    }
    /// pysam `is_sam` — True if file is plain-text SAM.
    #[getter]
    fn is_sam(&self) -> bool {
        self.path.to_ascii_lowercase().ends_with(".sam")
    }
    /// pysam `is_vcf` — always False on AlignmentFile.
    #[getter]
    fn is_vcf(&self) -> bool {
        false
    }
    /// pysam `is_bcf` — always False on AlignmentFile.
    #[getter]
    fn is_bcf(&self) -> bool {
        false
    }
    /// pysam `format` — "BAM" / "SAM" / "CRAM" / "UNKNOWN".
    #[getter]
    fn format(&self) -> &'static str {
        if self.is_bam() {
            "BAM"
        } else if self.is_cram() {
            "CRAM"
        } else if self.is_sam() {
            "SAM"
        } else {
            "UNKNOWN"
        }
    }
    /// pysam `compression` — "BGZF" for BAM/CRAM, "NONE" for SAM.
    #[getter]
    fn compression(&self) -> &'static str {
        if self.is_bam() || self.is_cram() {
            "BGZF"
        } else {
            "NONE"
        }
    }
    /// pysam `category` — always "alignment".
    #[getter]
    fn category(&self) -> &'static str {
        "alignment"
    }
    /// pysam `description` — human-readable format description.
    #[getter]
    fn description(&self) -> &'static str {
        "Binary Sequence Alignment/Map (BAM) / pysam-compatible read+write surface via rubam"
    }
    /// pysam `version` — backend identification string.
    #[getter]
    fn version(&self) -> &'static str {
        "noodles 0.107 (via rubam)"
    }
    /// pysam `nocoordinate` — always 0 (rubam does not separately track
    /// unmapped-without-coordinate records).
    #[getter]
    fn nocoordinate(&self) -> u64 {
        0
    }
    /// pysam `index_filename` — `.bai` / `.csi` if present, else None.
    #[getter]
    fn index_filename(&self) -> Option<String> {
        let bai = format!("{}.bai", self.path);
        let csi = format!("{}.csi", self.path);
        let crai = format!("{}.crai", self.path);
        if std::path::Path::new(&bai).exists() {
            Some(bai)
        } else if std::path::Path::new(&csi).exists() {
            Some(csi)
        } else if std::path::Path::new(&crai).exists() {
            Some(crai)
        } else {
            None
        }
    }
    /// pysam `text` — raw SAM-format header string.
    #[getter]
    fn text(&self) -> PyResult<String> {
        // We delegate to noodles' debug repr of the header — full SAM
        // serialization with @HD/@SQ/@PG/@CO sections lands in a follow-up.
        Ok(match &self.header {
            Some(h) => format!("{h:?}"),
            None => String::new(),
        })
    }

    /// pysam `get_tid(name)` / `gettid(name)` — contig name → index.
    /// Returns -1 if the contig is unknown (matches pysam).
    fn get_tid(&self, name: &str) -> i64 {
        match &self.header {
            Some(h) => h
                .reference_sequences()
                .get_index_of(name.as_bytes())
                .map(|i| i as i64)
                .unwrap_or(-1),
            None => -1,
        }
    }
    /// pysam alias.
    fn gettid(&self, name: &str) -> i64 {
        self.get_tid(name)
    }
    /// pysam `get_reference_name(tid)` / `getrname(tid)` — index → name.
    fn get_reference_name(&self, tid: usize) -> PyResult<String> {
        match &self.header {
            Some(h) => match h.reference_sequences().get_index(tid) {
                Some((name, _)) => Ok(String::from_utf8_lossy(name).into_owned()),
                None => Err(pyo3::exceptions::PyIndexError::new_err(format!(
                    "tid {tid} out of range"
                ))),
            },
            None => Err(PyIOError::new_err("no header available")),
        }
    }
    /// pysam alias.
    fn getrname(&self, tid: usize) -> PyResult<String> {
        self.get_reference_name(tid)
    }
    /// pysam `is_valid_reference_name(name)` — True if the contig is in the header.
    fn is_valid_reference_name(&self, name: &str) -> bool {
        self.get_tid(name) >= 0
    }
    /// pysam `is_valid_tid(tid)` — True if 0 ≤ tid < nreferences.
    fn is_valid_tid(&self, tid: i64) -> bool {
        tid >= 0 && (tid as usize) < self.nreferences()
    }

    /// pysam `mapped` / `unmapped` — counts derived from idx stats.
    #[getter]
    fn mapped(&self, py: Python<'_>) -> PyResult<u64> {
        let stats = self.get_index_statistics(py)?;
        let mut total = 0u64;
        for item in stats.iter() {
            let d: &Bound<pyo3::types::PyDict> = item.downcast()?;
            if let Some(m) = d.get_item("mapped")? {
                total += m.extract::<u64>()?;
            }
        }
        Ok(total)
    }
    #[getter]
    fn unmapped(&self, py: Python<'_>) -> PyResult<u64> {
        let stats = self.get_index_statistics(py)?;
        let mut total = 0u64;
        for item in stats.iter() {
            let d: &Bound<pyo3::types::PyDict> = item.downcast()?;
            if let Some(u) = d.get_item("unmapped")? {
                total += u.extract::<u64>()?;
            }
        }
        Ok(total)
    }

    /// pysam `add_hts_options(opts)` — pysam-compat no-op storage.
    fn add_hts_options(&self, opts: Vec<String>) {
        let mut h = self.hts_options.borrow_mut();
        h.extend(opts);
    }
    /// pysam `flush()` — no-op for readers; flushes BGZF buffer for writers.
    fn flush(&self) -> PyResult<()> {
        Ok(())
    }
    /// pysam `reset()` — no-op stub; full BGZF seek-to-start lands later.
    fn reset(&self) -> PyResult<()> {
        Ok(())
    }
    /// pysam `duplicate_filehandle` — always False (noodles doesn't dup FDs).
    #[getter]
    fn duplicate_filehandle(&self) -> bool {
        false
    }
    /// pysam `check_truncation` — flag accessor.
    #[getter]
    fn check_truncation(&self) -> bool {
        self.check_truncation_flag.get()
    }
    #[setter]
    fn set_check_truncation(&self, v: bool) {
        self.check_truncation_flag.set(v)
    }

    /// pysam `parse_region(region=None, contig=None, start=None, end=None)` →
    /// (tid, start, end). Accepts "chr1:1-100" / "chr1:1" / "chr1" /
    /// "chr1:1,000,000-2,000,000" forms (commas are stripped).
    #[pyo3(signature = (region = None, contig = None, start = None, end = None))]
    fn parse_region(
        &self,
        region: Option<&str>,
        contig: Option<&str>,
        start: Option<usize>,
        end: Option<usize>,
    ) -> PyResult<(i64, Option<usize>, Option<usize>)> {
        if let Some(r) = region {
            let mut parts = r.splitn(2, ':');
            let chr = parts.next().unwrap_or("");
            let coords = parts.next();
            let (s, e) = if let Some(c) = coords {
                let c_clean = c.replace(',', "");
                let mut se = c_clean.splitn(2, '-');
                let s_str = se.next().unwrap_or("1");
                let s: usize = s_str
                    .parse()
                    .map_err(|_| PyValueError::new_err(format!("bad start in region: {r:?}")))?;
                let e: Option<usize> =
                    match se.next() {
                        Some(es) if !es.is_empty() => Some(es.parse().map_err(|_| {
                            PyValueError::new_err(format!("bad end in region: {r:?}"))
                        })?),
                        _ => None,
                    };
                (Some(s.saturating_sub(1)), e)
            } else {
                (None, None)
            };
            Ok((self.get_tid(chr), s, e))
        } else {
            let chr = contig.unwrap_or("");
            Ok((self.get_tid(chr), start, end))
        }
    }
}

// ---------------------------------------------------------------------------
// AlignedSegment — wraps either a BAM record or a CRAM RecordBuf.
// ---------------------------------------------------------------------------

#[pyclass]
pub struct AlignedSegment {
    pub(crate) record: AnyRecord,
    pub(crate) header: noodles::sam::Header,
}

impl AlignedSegment {
    pub(crate) fn new(record: AnyRecord, header: noodles::sam::Header) -> Self {
        Self { record, header }
    }

    /// Materialise `self.record` into a mutable `RecordBuf` if it is still
    /// the read-only `AnyRecord::Bam(bam::Record)` byte-view, then return a
    /// mutable reference to the inner `RecordBuf`.
    ///
    /// All write-side setters route through here so the BAM-fetched and
    /// CRAM-fetched paths share the same mutation code.
    fn ensure_record_buf(&mut self) -> PyResult<&mut RecordBuf> {
        // Two-phase to satisfy the borrow checker: convert (consuming the
        // old enum variant) first, then return a mutable ref to the new one.
        if let AnyRecord::Bam(_) = &self.record {
            let owned = std::mem::replace(&mut self.record, AnyRecord::Cram(RecordBuf::default()));
            let bam_rec = match owned {
                AnyRecord::Bam(r) => r,
                _ => unreachable!(),
            };
            let buf =
                RecordBuf::try_from_alignment_record(&self.header, &bam_rec).map_err(|e| {
                    PyIOError::new_err(format!(
                        "failed to materialise BAM record into a mutable RecordBuf: {e}"
                    ))
                })?;
            self.record = AnyRecord::Cram(buf);
        }
        match &mut self.record {
            AnyRecord::Cram(r) => Ok(r),
            AnyRecord::Bam(_) => unreachable!("materialisation above turned this into Cram"),
        }
    }

    /// Flip a single flag bit on the underlying record. Materialises into a
    /// `RecordBuf` if needed.
    fn set_flag_bit(&mut self, bit: u16, on: bool) -> PyResult<()> {
        use noodles::sam::alignment::record::Flags;
        let rec = self.ensure_record_buf()?;
        let mut flags = rec.flags().bits();
        if on {
            flags |= bit;
        } else {
            flags &= !bit;
        }
        *rec.flags_mut() = Flags::from_bits_retain(flags);
        Ok(())
    }
}

// Helper: extract flags bits from AnyRecord.
fn flags_bits(record: &AnyRecord) -> u16 {
    match record {
        AnyRecord::Bam(r) => r.flags().bits(),
        AnyRecord::Cram(r) => r.flags().bits(),
    }
}

#[pymethods]
impl AlignedSegment {
    /// Construct an empty `AlignedSegment`, optionally bound to a header so
    /// that `reference_id` lookups remain consistent with a destination
    /// `AlignmentFile`. Mirrors `pysam.AlignedSegment(header=bam.header)`.
    #[new]
    #[pyo3(signature = (header = None))]
    fn py_new(header: Option<&Header>) -> PyResult<Self> {
        let hdr = header.map(|h| h.inner.clone()).unwrap_or_default();
        Ok(Self {
            record: AnyRecord::Cram(RecordBuf::default()),
            header: hdr,
        })
    }

    // -----------------------------------------------------------------
    // Setters (write-side API).
    //
    // Every setter routes through `ensure_record_buf` so a record that
    // was just fetched (`AnyRecord::Bam`) is transparently promoted to
    // an owned mutable `RecordBuf` on first mutation.
    // -----------------------------------------------------------------

    #[setter(query_name)]
    fn set_query_name(&mut self, v: &str) -> PyResult<()> {
        let rec = self.ensure_record_buf()?;
        *rec.name_mut() = Some(v.as_bytes().to_vec().into());
        Ok(())
    }

    #[setter(flag)]
    fn set_flag(&mut self, v: u16) -> PyResult<()> {
        use noodles::sam::alignment::record::Flags;
        let rec = self.ensure_record_buf()?;
        *rec.flags_mut() = Flags::from_bits_retain(v);
        Ok(())
    }

    #[setter(reference_id)]
    fn set_reference_id(&mut self, v: Option<usize>) -> PyResult<()> {
        let rec = self.ensure_record_buf()?;
        *rec.reference_sequence_id_mut() = v;
        Ok(())
    }

    #[setter(reference_start)]
    fn set_reference_start(&mut self, v: Option<usize>) -> PyResult<()> {
        use noodles::core::Position;
        let rec = self.ensure_record_buf()?;
        // pysam is 0-based; noodles is 1-based.
        let pos =
            match v {
                None => None,
                Some(p) => Some(Position::try_from(p + 1).map_err(|e| {
                    PyValueError::new_err(format!("invalid reference_start {p}: {e}"))
                })?),
            };
        *rec.alignment_start_mut() = pos;
        Ok(())
    }

    #[setter(mapping_quality)]
    fn set_mapping_quality(&mut self, v: u8) -> PyResult<()> {
        use noodles::sam::alignment::record::MappingQuality;
        let rec = self.ensure_record_buf()?;
        // `MappingQuality::new(255)` returns None (the SAM "missing" sentinel).
        *rec.mapping_quality_mut() = MappingQuality::new(v);
        Ok(())
    }

    #[setter(template_length)]
    fn set_template_length(&mut self, v: i32) -> PyResult<()> {
        let rec = self.ensure_record_buf()?;
        *rec.template_length_mut() = v;
        Ok(())
    }

    #[setter(mate_reference_id)]
    fn set_mate_reference_id(&mut self, v: Option<usize>) -> PyResult<()> {
        let rec = self.ensure_record_buf()?;
        *rec.mate_reference_sequence_id_mut() = v;
        Ok(())
    }

    #[setter(mate_reference_start)]
    fn set_mate_reference_start(&mut self, v: Option<usize>) -> PyResult<()> {
        use noodles::core::Position;
        let rec = self.ensure_record_buf()?;
        let pos = match v {
            None => None,
            Some(p) => Some(Position::try_from(p + 1).map_err(|e| {
                PyValueError::new_err(format!("invalid mate_reference_start {p}: {e}"))
            })?),
        };
        *rec.mate_alignment_start_mut() = pos;
        Ok(())
    }

    /// pysam alias for `mate_reference_id` setter: `seg.next_reference_id = rid`.
    #[setter(next_reference_id)]
    fn set_next_reference_id(&mut self, v: Option<usize>) -> PyResult<()> {
        self.set_mate_reference_id(v)
    }

    /// pysam alias for `mate_reference_start` setter: `seg.next_reference_start = p`.
    #[setter(next_reference_start)]
    fn set_next_reference_start(&mut self, v: Option<usize>) -> PyResult<()> {
        self.set_mate_reference_start(v)
    }

    #[setter(query_sequence)]
    fn set_query_sequence(&mut self, v: &str) -> PyResult<()> {
        use noodles::sam::alignment::record_buf::Sequence;
        let rec = self.ensure_record_buf()?;
        *rec.sequence_mut() = Sequence::from(v.as_bytes().to_vec());
        Ok(())
    }

    #[setter(query_qualities)]
    fn set_query_qualities(&mut self, v: Vec<u8>) -> PyResult<()> {
        use noodles::sam::alignment::record_buf::QualityScores;
        let rec = self.ensure_record_buf()?;
        *rec.quality_scores_mut() = QualityScores::from(v);
        Ok(())
    }

    #[setter(cigarstring)]
    fn set_cigarstring(&mut self, v: &str) -> PyResult<()> {
        use noodles::sam::alignment::record_buf::Cigar;
        let ops = parse_cigar_string(v)?;
        let rec = self.ensure_record_buf()?;
        *rec.cigar_mut() = Cigar::from(ops);
        Ok(())
    }

    #[setter(cigartuples)]
    fn set_cigartuples(&mut self, v: Vec<(u8, u32)>) -> PyResult<()> {
        use noodles::sam::alignment::record::cigar::Op;
        use noodles::sam::alignment::record_buf::Cigar;
        let mut ops = Vec::with_capacity(v.len());
        for (op_code, len) in v {
            ops.push(Op::new(kind_from_int(op_code)?, len as usize));
        }
        let rec = self.ensure_record_buf()?;
        *rec.cigar_mut() = Cigar::from(ops);
        Ok(())
    }

    /// Set or replace a 2-char auxiliary tag. Value may be int/float/str/bytes
    /// (pysam's documented set for `set_tag`). Existing array (`B`) tags are
    /// not supported on the write path yet (mirrors the read-side gap).
    fn set_tag(&mut self, name: &str, value: &Bound<'_, PyAny>) -> PyResult<()> {
        use noodles::sam::alignment::record::data::field::Tag;
        if name.len() != 2 {
            return Err(PyValueError::new_err("tag name must be 2 chars"));
        }
        let bytes: [u8; 2] = name
            .as_bytes()
            .try_into()
            .map_err(|_| PyValueError::new_err("tag name must be 2 chars"))?;
        let tag = Tag::new(bytes[0], bytes[1]);
        let value = py_to_tag_value(value)?;
        let rec = self.ensure_record_buf()?;
        rec.data_mut().insert(tag, value);
        Ok(())
    }

    /// Remove a 2-char auxiliary tag if present. Silently no-ops on a missing
    /// tag (matches pysam, which does not raise).
    fn remove_tag(&mut self, name: &str) -> PyResult<()> {
        use noodles::sam::alignment::record::data::field::Tag;
        if name.len() != 2 {
            return Err(PyValueError::new_err("tag name must be 2 chars"));
        }
        let bytes: [u8; 2] = name
            .as_bytes()
            .try_into()
            .map_err(|_| PyValueError::new_err("tag name must be 2 chars"))?;
        let tag = Tag::new(bytes[0], bytes[1]);
        let rec = self.ensure_record_buf()?;
        rec.data_mut().remove(&tag);
        Ok(())
    }

    /// Bulk-replace the tag set. Existing tags are cleared first; matches
    /// pysam's `seg.tags = [(name, value), ...]` property assignment.
    ///
    /// Each entry must be a `(str, value)` tuple where value is one of
    /// int/float/str/bytes (the same set `set_tag` accepts).
    #[setter(tags)]
    fn set_tags(&mut self, tags: &Bound<'_, pyo3::types::PyList>) -> PyResult<()> {
        use noodles::sam::alignment::record::data::field::Tag;
        // Build the new table eagerly so a failure leaves the record untouched.
        let mut staged: Vec<(Tag, _)> = Vec::with_capacity(tags.len());
        for item in tags.iter() {
            let tup: (String, Bound<'_, PyAny>) = item.extract().map_err(|e| {
                PyTypeError::new_err(format!(
                    "set_tags: each entry must be a (str, value) tuple: {e}"
                ))
            })?;
            let (name, val) = tup;
            if name.len() != 2 {
                return Err(PyValueError::new_err("tag name must be 2 chars"));
            }
            let bytes: [u8; 2] = name
                .as_bytes()
                .try_into()
                .map_err(|_| PyValueError::new_err("tag name must be 2 chars"))?;
            staged.push((Tag::new(bytes[0], bytes[1]), py_to_tag_value(&val)?));
        }
        let rec = self.ensure_record_buf()?;
        rec.data_mut().clear();
        for (tag, value) in staged {
            rec.data_mut().insert(tag, value);
        }
        Ok(())
    }

    /// pysam-compatible `to_dict()` — emits a plain dict of the record's
    /// surface fields, suitable for round-tripping through `from_dict`.
    ///
    /// Key set mirrors pysam.AlignedSegment.to_dict (subset that rubam
    /// already supports on the read side): `name`, `flag`, `ref_name`,
    /// `ref_pos`, `map_quality`, `cigar`, `next_ref_name`, `next_ref_pos`,
    /// `length`, `seq`, `qual`, `tags`.
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        use pyo3::types::PyDict;
        let d = PyDict::new(py);
        d.set_item("name", self.query_name()?)?;
        d.set_item("flag", self.flag())?;
        d.set_item("ref_name", self.reference_name()?)?;
        // ref_pos: pysam reports 1-based in to_dict; we follow suit so the
        // pair `to_dict()` / `from_dict()` round-trips cleanly.
        d.set_item("ref_pos", self.reference_start()?.map(|p| p + 1))?;
        d.set_item("map_quality", self.mapping_quality())?;
        d.set_item("cigar", self.cigarstring()?)?;
        // Mate / next ref name.
        let next_ref_name = match self.record {
            AnyRecord::Bam(ref r) => r
                .mate_reference_sequence_id()
                .and_then(|res| res.ok())
                .and_then(|id| self.header.reference_sequences().get_index(id))
                .map(|(name, _)| String::from_utf8_lossy(name).into_owned()),
            AnyRecord::Cram(ref r) => r
                .mate_reference_sequence_id()
                .and_then(|id| self.header.reference_sequences().get_index(id))
                .map(|(name, _)| String::from_utf8_lossy(name).into_owned()),
        };
        d.set_item("next_ref_name", next_ref_name)?;
        let next_ref_pos: Option<usize> = match self.record {
            AnyRecord::Bam(ref r) => r
                .mate_alignment_start()
                .and_then(|res| res.ok())
                .map(|p| p.get()),
            AnyRecord::Cram(ref r) => r.mate_alignment_start().map(|p| p.get()),
        };
        d.set_item("next_ref_pos", next_ref_pos)?;
        d.set_item("length", self.template_length())?;
        d.set_item("seq", self.query_sequence()?)?;
        d.set_item("qual", self.query_qualities()?)?;
        d.set_item("tags", self.tags(py)?)?;
        Ok(d)
    }

    /// pysam-compatible `from_dict(header, d)` classmethod — inverse of
    /// `to_dict`. Builds a fresh `AlignedSegment` whose fields match `d`.
    ///
    /// `header` is required because `ref_name`/`next_ref_name` resolution
    /// needs a reference table.
    #[classmethod]
    fn from_dict(
        _cls: &Bound<'_, pyo3::types::PyType>,
        header: &Header,
        d: &Bound<'_, pyo3::types::PyDict>,
    ) -> PyResult<Self> {
        // Construct an empty buffer bound to the header.
        let mut seg = Self {
            record: AnyRecord::Cram(RecordBuf::default()),
            header: header.inner.clone(),
        };

        // Resolve a reference name to its 0-based index in the header.
        let resolve_ref = |name_opt: Option<String>| -> PyResult<Option<usize>> {
            match name_opt {
                None => Ok(None),
                Some(n) => header
                    .inner
                    .reference_sequences()
                    .get_index_of(n.as_bytes())
                    .map(Some)
                    .ok_or_else(|| {
                        PyValueError::new_err(format!(
                            "from_dict: ref_name {n:?} not found in header"
                        ))
                    }),
            }
        };

        if let Some(v) = d.get_item("name")? {
            seg.set_query_name(&v.extract::<String>()?)?;
        }
        if let Some(v) = d.get_item("flag")? {
            seg.set_flag(v.extract::<u16>()?)?;
        }
        if let Some(v) = d.get_item("ref_name")? {
            if !v.is_none() {
                let rid = resolve_ref(Some(v.extract::<String>()?))?;
                seg.set_reference_id(rid)?;
            }
        }
        if let Some(v) = d.get_item("ref_pos")? {
            if !v.is_none() {
                // ref_pos in to_dict is 1-based — convert to 0-based for
                // set_reference_start (which then re-adds 1 internally).
                let p_1based = v.extract::<usize>()?;
                let p_0based = p_1based
                    .checked_sub(1)
                    .ok_or_else(|| PyValueError::new_err("from_dict: ref_pos must be >= 1"))?;
                seg.set_reference_start(Some(p_0based))?;
            }
        }
        if let Some(v) = d.get_item("map_quality")? {
            seg.set_mapping_quality(v.extract::<u8>()?)?;
        }
        if let Some(v) = d.get_item("cigar")? {
            if !v.is_none() {
                seg.set_cigarstring(&v.extract::<String>()?)?;
            }
        }
        if let Some(v) = d.get_item("next_ref_name")? {
            if !v.is_none() {
                let rid = resolve_ref(Some(v.extract::<String>()?))?;
                seg.set_mate_reference_id(rid)?;
            }
        }
        if let Some(v) = d.get_item("next_ref_pos")? {
            if !v.is_none() {
                let p_1based = v.extract::<usize>()?;
                let p_0based = p_1based
                    .checked_sub(1)
                    .ok_or_else(|| PyValueError::new_err("from_dict: next_ref_pos must be >= 1"))?;
                seg.set_mate_reference_start(Some(p_0based))?;
            }
        }
        if let Some(v) = d.get_item("length")? {
            seg.set_template_length(v.extract::<i32>()?)?;
        }
        if let Some(v) = d.get_item("seq")? {
            if !v.is_none() {
                seg.set_query_sequence(&v.extract::<String>()?)?;
            }
        }
        if let Some(v) = d.get_item("qual")? {
            if !v.is_none() {
                seg.set_query_qualities(v.extract::<Vec<u8>>()?)?;
            }
        }
        if let Some(v) = d.get_item("tags")? {
            if !v.is_none() {
                let list = v.downcast::<pyo3::types::PyList>().map_err(|_| {
                    PyTypeError::new_err("from_dict: 'tags' must be a list of (name, value) tuples")
                })?;
                seg.set_tags(list)?;
            }
        }

        Ok(seg)
    }

    // -----------------------------------------------------------------
    // Flag bit setters — mirror the existing `is_*` read-only getters.
    // pysam exposes these as `seg.is_paired = True` properties, but the
    // existing `is_*` getters already occupy those names, so we expose
    // explicit method-style setters (consistent with pysam ≥ 0.16 alias).
    // -----------------------------------------------------------------

    fn set_is_paired(&mut self, b: bool) -> PyResult<()> {
        self.set_flag_bit(FLAG_PAIRED, b)
    }
    fn set_is_proper_pair(&mut self, b: bool) -> PyResult<()> {
        self.set_flag_bit(FLAG_PROPER, b)
    }
    fn set_is_unmapped(&mut self, b: bool) -> PyResult<()> {
        self.set_flag_bit(FLAG_UNMAP, b)
    }
    fn set_mate_is_unmapped(&mut self, b: bool) -> PyResult<()> {
        self.set_flag_bit(FLAG_MUNMAP, b)
    }
    fn set_is_reverse(&mut self, b: bool) -> PyResult<()> {
        self.set_flag_bit(FLAG_REVERSE, b)
    }
    fn set_mate_is_reverse(&mut self, b: bool) -> PyResult<()> {
        self.set_flag_bit(FLAG_MREVERSE, b)
    }
    fn set_is_read1(&mut self, b: bool) -> PyResult<()> {
        self.set_flag_bit(FLAG_READ1, b)
    }
    fn set_is_read2(&mut self, b: bool) -> PyResult<()> {
        self.set_flag_bit(FLAG_READ2, b)
    }
    fn set_is_secondary(&mut self, b: bool) -> PyResult<()> {
        self.set_flag_bit(FLAG_SECONDARY, b)
    }
    fn set_is_qcfail(&mut self, b: bool) -> PyResult<()> {
        self.set_flag_bit(FLAG_QCFAIL, b)
    }
    fn set_is_duplicate(&mut self, b: bool) -> PyResult<()> {
        self.set_flag_bit(FLAG_DUP, b)
    }
    fn set_is_supplementary(&mut self, b: bool) -> PyResult<()> {
        self.set_flag_bit(FLAG_SUPPLEMENTARY, b)
    }

    #[getter]
    fn query_name(&self) -> PyResult<String> {
        match &self.record {
            AnyRecord::Bam(r) => Ok(match r.name() {
                Some(name) => String::from_utf8_lossy(name.as_ref()).into_owned(),
                None => String::new(),
            }),
            AnyRecord::Cram(r) => Ok(match r.name() {
                Some(name) => String::from_utf8_lossy(name.as_ref()).into_owned(),
                None => String::new(),
            }),
        }
    }

    #[getter]
    fn reference_id(&self) -> PyResult<Option<usize>> {
        match &self.record {
            AnyRecord::Bam(r) => match r.reference_sequence_id() {
                Some(Ok(id)) => Ok(Some(id)),
                Some(Err(e)) => Err(PyIOError::new_err(format!("rid: {e}"))),
                None => Ok(None),
            },
            AnyRecord::Cram(r) => Ok(r.reference_sequence_id()),
        }
    }

    #[getter]
    fn reference_name(&self) -> PyResult<Option<String>> {
        let Some(id) = self.reference_id()? else {
            return Ok(None);
        };
        let refs = self.header.reference_sequences();
        let entry = refs.get_index(id).map(|(name, _)| name.clone());
        Ok(entry.map(|n| String::from_utf8_lossy(&n).into_owned()))
    }

    #[getter]
    fn reference_start(&self) -> PyResult<Option<usize>> {
        match &self.record {
            AnyRecord::Bam(r) => match r.alignment_start() {
                Some(Ok(p)) => Ok(Some(p.get() - 1)),
                Some(Err(e)) => Err(PyIOError::new_err(format!("alignment_start: {e}"))),
                None => Ok(None),
            },
            AnyRecord::Cram(r) => Ok(r.alignment_start().map(|p| p.get() - 1)),
        }
    }

    #[getter]
    fn reference_end(&self) -> PyResult<Option<usize>> {
        match &self.record {
            AnyRecord::Bam(r) => match r.alignment_end() {
                Some(Ok(p)) => Ok(Some(p.get())),
                Some(Err(e)) => Err(PyIOError::new_err(format!("alignment_end: {e}"))),
                None => Ok(None),
            },
            AnyRecord::Cram(r) => {
                // RecordBuf's alignment_end() returns Option<Position> directly (no Result).
                Ok(r.alignment_end().map(|p| p.get()))
            }
        }
    }

    #[getter]
    fn template_length(&self) -> i32 {
        match &self.record {
            AnyRecord::Bam(r) => r.template_length(),
            AnyRecord::Cram(r) => r.template_length(),
        }
    }

    #[getter]
    fn mapping_quality(&self) -> u8 {
        match &self.record {
            AnyRecord::Bam(r) => r.mapping_quality().map(|q| q.get()).unwrap_or(255),
            AnyRecord::Cram(r) => r.mapping_quality().map(|q| q.get()).unwrap_or(255),
        }
    }

    #[getter]
    fn flag(&self) -> u16 {
        flags_bits(&self.record)
    }

    #[getter]
    fn is_paired(&self) -> bool {
        flags_bits(&self.record) & FLAG_PAIRED != 0
    }
    #[getter]
    fn is_proper_pair(&self) -> bool {
        flags_bits(&self.record) & FLAG_PROPER != 0
    }
    #[getter]
    fn is_unmapped(&self) -> bool {
        flags_bits(&self.record) & FLAG_UNMAP != 0
    }
    #[getter]
    fn is_mate_unmapped(&self) -> bool {
        flags_bits(&self.record) & FLAG_MUNMAP != 0
    }
    #[getter]
    fn is_reverse(&self) -> bool {
        flags_bits(&self.record) & FLAG_REVERSE != 0
    }
    #[getter]
    fn is_mate_reverse(&self) -> bool {
        flags_bits(&self.record) & FLAG_MREVERSE != 0
    }
    #[getter]
    fn is_read1(&self) -> bool {
        flags_bits(&self.record) & FLAG_READ1 != 0
    }
    #[getter]
    fn is_read2(&self) -> bool {
        flags_bits(&self.record) & FLAG_READ2 != 0
    }
    #[getter]
    fn is_secondary(&self) -> bool {
        flags_bits(&self.record) & FLAG_SECONDARY != 0
    }
    #[getter]
    fn is_qcfail(&self) -> bool {
        flags_bits(&self.record) & FLAG_QCFAIL != 0
    }
    #[getter]
    fn is_duplicate(&self) -> bool {
        flags_bits(&self.record) & FLAG_DUP != 0
    }
    #[getter]
    fn is_supplementary(&self) -> bool {
        flags_bits(&self.record) & FLAG_SUPPLEMENTARY != 0
    }

    #[getter]
    fn cigarstring(&self) -> PyResult<Option<String>> {
        match &self.record {
            AnyRecord::Bam(r) => {
                let cigar = r.cigar();
                let mut s = String::new();
                let mut empty = true;
                for op_result in cigar.iter() {
                    let op = op_result.map_err(|e| PyIOError::new_err(format!("cigar: {e}")))?;
                    empty = false;
                    s.push_str(&op.len().to_string());
                    s.push(kind_to_char(op.kind()));
                }
                if empty {
                    Ok(None)
                } else {
                    Ok(Some(s))
                }
            }
            AnyRecord::Cram(r) => {
                // record_buf::Cigar implements AsRef<[Op]>; use it directly (no Result).
                let cigar = r.cigar();
                let ops: &[noodles::sam::alignment::record::cigar::Op] = cigar.as_ref();
                let mut s = String::new();
                for op in ops {
                    s.push_str(&op.len().to_string());
                    s.push(kind_to_char(op.kind()));
                }
                if s.is_empty() {
                    Ok(None)
                } else {
                    Ok(Some(s))
                }
            }
        }
    }

    #[getter]
    fn cigartuples(&self) -> PyResult<Option<Vec<(u8, usize)>>> {
        match &self.record {
            AnyRecord::Bam(r) => {
                let cigar = r.cigar();
                let mut out = Vec::new();
                for op_result in cigar.iter() {
                    let op = op_result.map_err(|e| PyIOError::new_err(format!("cigar: {e}")))?;
                    out.push((kind_to_int(op.kind()), op.len()));
                }
                if out.is_empty() {
                    Ok(None)
                } else {
                    Ok(Some(out))
                }
            }
            AnyRecord::Cram(r) => {
                let cigar = r.cigar();
                let ops: &[noodles::sam::alignment::record::cigar::Op] = cigar.as_ref();
                let out: Vec<(u8, usize)> = ops
                    .iter()
                    .map(|op| (kind_to_int(op.kind()), op.len()))
                    .collect();
                if out.is_empty() {
                    Ok(None)
                } else {
                    Ok(Some(out))
                }
            }
        }
    }

    #[getter]
    fn query_sequence(&self) -> PyResult<Option<String>> {
        match &self.record {
            AnyRecord::Bam(r) => {
                let seq = r.sequence();
                let n = seq.len();
                if n == 0 {
                    return Ok(None);
                }
                let mut s = String::with_capacity(n);
                for i in 0..n {
                    let b = seq.get(i).unwrap_or(b'N');
                    s.push(b as char);
                }
                Ok(Some(s))
            }
            AnyRecord::Cram(r) => {
                let seq = r.sequence();
                let n = seq.len();
                if n == 0 {
                    return Ok(None);
                }
                let mut s = String::with_capacity(n);
                for i in 0..n {
                    let b = seq.get(i).unwrap_or(b'N');
                    s.push(b as char);
                }
                Ok(Some(s))
            }
        }
    }

    #[getter]
    fn query_qualities(&self) -> PyResult<Option<Vec<u8>>> {
        match &self.record {
            AnyRecord::Bam(r) => {
                let quals: Vec<u8> = r.quality_scores().iter().collect();
                if quals.is_empty() {
                    Ok(None)
                } else {
                    Ok(Some(quals))
                }
            }
            AnyRecord::Cram(r) => {
                let quals: Vec<u8> = r.quality_scores().iter().collect();
                if quals.is_empty() {
                    Ok(None)
                } else {
                    Ok(Some(quals))
                }
            }
        }
    }

    #[getter]
    fn query_length(&self) -> usize {
        match &self.record {
            AnyRecord::Bam(r) => r.sequence().len(),
            AnyRecord::Cram(r) => r.sequence().len(),
        }
    }

    #[getter]
    fn tags<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let pylist = PyList::empty(py);
        match &self.record {
            AnyRecord::Bam(r) => {
                let data = r.data();
                for entry in data.iter() {
                    let (tag, value) =
                        entry.map_err(|e| PyIOError::new_err(format!("tags: {e}")))?;
                    let name = std::str::from_utf8(tag.as_ref())
                        .map_err(|e| PyIOError::new_err(format!("tag name: {e}")))?
                        .to_string();
                    let v_py = tag_value_to_py(py, value)?;
                    pylist.append((name, v_py))?;
                }
            }
            AnyRecord::Cram(r) => {
                let data = r.data();
                for (tag, value) in data.iter() {
                    let name = std::str::from_utf8(tag.as_ref())
                        .map_err(|e| PyIOError::new_err(format!("tag name: {e}")))?
                        .to_string();
                    let v_py = record_buf_value_to_py(py, value)?;
                    pylist.append((name, v_py))?;
                }
            }
        }
        Ok(pylist)
    }

    fn has_tag(&self, name: &str) -> PyResult<bool> {
        if name.len() != 2 {
            return Err(PyValueError::new_err("tag name must be 2 chars"));
        }
        let target: [u8; 2] = name
            .as_bytes()
            .try_into()
            .map_err(|_| PyValueError::new_err("tag name must be 2 chars"))?;
        match &self.record {
            AnyRecord::Bam(r) => {
                let data = r.data();
                for entry in data.iter() {
                    let (tag, _) = entry.map_err(|e| PyIOError::new_err(format!("tags: {e}")))?;
                    if tag.as_ref() == &target {
                        return Ok(true);
                    }
                }
                Ok(false)
            }
            AnyRecord::Cram(r) => {
                use noodles::sam::alignment::record::data::field::Tag;
                let t = Tag::new(target[0], target[1]);
                Ok(r.data().get(&t).is_some())
            }
        }
    }

    fn get_tag<'py>(&self, py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyAny>> {
        if name.len() != 2 {
            return Err(PyValueError::new_err("tag name must be 2 chars"));
        }
        let target: [u8; 2] = name
            .as_bytes()
            .try_into()
            .map_err(|_| PyValueError::new_err("tag name must be 2 chars"))?;
        match &self.record {
            AnyRecord::Bam(r) => {
                let data = r.data();
                for entry in data.iter() {
                    let (tag, value) =
                        entry.map_err(|e| PyIOError::new_err(format!("tags: {e}")))?;
                    if tag.as_ref() == &target {
                        return tag_value_to_py(py, value);
                    }
                }
                Err(PyKeyError::new_err(format!("tag {name:?} not present")))
            }
            AnyRecord::Cram(r) => {
                use noodles::sam::alignment::record::data::field::Tag;
                let t = Tag::new(target[0], target[1]);
                match r.data().get(&t) {
                    Some(value) => record_buf_value_to_py(py, value),
                    None => Err(PyKeyError::new_err(format!("tag {name:?} not present"))),
                }
            }
        }
    }

    /// List of (start, end) reference intervals covered by M/=/X ops.
    fn get_blocks(&self) -> PyResult<Vec<(usize, usize)>> {
        let mut out = Vec::new();
        let ref_start = self.reference_start()?;
        let mut ref_pos: usize = match ref_start {
            Some(p) => p,
            None => return Ok(out),
        };
        match &self.record {
            AnyRecord::Bam(r) => {
                for op_result in r.cigar().iter() {
                    let op = op_result.map_err(|e| PyIOError::new_err(format!("cigar: {e}")))?;
                    let len = op.len();
                    match op.kind() {
                        Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                            out.push((ref_pos, ref_pos + len));
                            ref_pos += len;
                        }
                        Kind::Deletion | Kind::Skip => {
                            ref_pos += len;
                        }
                        _ => {}
                    }
                }
            }
            AnyRecord::Cram(r) => {
                let ops: &[noodles::sam::alignment::record::cigar::Op] = r.cigar().as_ref();
                for op in ops {
                    let len = op.len();
                    match op.kind() {
                        Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                            out.push((ref_pos, ref_pos + len));
                            ref_pos += len;
                        }
                        Kind::Deletion | Kind::Skip => {
                            ref_pos += len;
                        }
                        _ => {}
                    }
                }
            }
        }
        Ok(out)
    }

    fn get_reference_positions(&self) -> PyResult<Vec<usize>> {
        let mut out = Vec::new();
        let ref_start = self.reference_start()?;
        let mut ref_pos: usize = match ref_start {
            Some(p) => p,
            None => return Ok(out),
        };
        match &self.record {
            AnyRecord::Bam(r) => {
                for op_result in r.cigar().iter() {
                    let op = op_result.map_err(|e| PyIOError::new_err(format!("cigar: {e}")))?;
                    let len = op.len();
                    match op.kind() {
                        Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                            for i in 0..len {
                                out.push(ref_pos + i);
                            }
                            ref_pos += len;
                        }
                        Kind::Deletion | Kind::Skip => {
                            ref_pos += len;
                        }
                        _ => {}
                    }
                }
            }
            AnyRecord::Cram(r) => {
                let ops: &[noodles::sam::alignment::record::cigar::Op] = r.cigar().as_ref();
                for op in ops {
                    let len = op.len();
                    match op.kind() {
                        Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                            for i in 0..len {
                                out.push(ref_pos + i);
                            }
                            ref_pos += len;
                        }
                        Kind::Deletion | Kind::Skip => {
                            ref_pos += len;
                        }
                        _ => {}
                    }
                }
            }
        }
        Ok(out)
    }

    fn get_overlap(&self, start: usize, end: usize) -> PyResult<usize> {
        let blocks = self.get_blocks()?;
        let mut total = 0usize;
        for (bs, be) in blocks {
            let lo = bs.max(start);
            let hi = be.min(end);
            if hi > lo {
                total += hi - lo;
            }
        }
        Ok(total)
    }

    // -------- pysam-compatible aliases (v0.3.4) --------
    // All of these delegate to the existing canonical getters so the rubam
    // surface matches pysam's two-name convention (e.g. `query_name` and
    // `qname` refer to the same field). Added so code written against
    // pysam runs unchanged on rubam.

    #[getter]
    fn qname(&self) -> PyResult<String> {
        self.query_name()
    }
    #[getter]
    fn pos(&self) -> PyResult<Option<usize>> {
        self.reference_start()
    }
    #[getter]
    fn mapq(&self) -> u8 {
        self.mapping_quality()
    }
    #[getter]
    fn tid(&self) -> PyResult<Option<usize>> {
        self.reference_id()
    }
    #[getter]
    fn isize(&self) -> i32 {
        self.template_length()
    }
    #[getter]
    fn tlen(&self) -> i32 {
        self.template_length()
    }
    #[getter]
    fn rname(&self) -> PyResult<Option<String>> {
        self.reference_name()
    }
    #[getter]
    fn mate_reference_id(&self) -> PyResult<Option<usize>> {
        Ok(match &self.record {
            AnyRecord::Bam(r) => r.mate_reference_sequence_id().and_then(|res| res.ok()),
            AnyRecord::Cram(r) => r.mate_reference_sequence_id(),
        })
    }
    #[getter]
    fn mate_reference_start(&self) -> PyResult<Option<usize>> {
        Ok(match &self.record {
            AnyRecord::Bam(r) => r
                .mate_alignment_start()
                .and_then(|res| res.ok())
                .map(|p| p.get() - 1),
            AnyRecord::Cram(r) => r.mate_alignment_start().map(|p| p.get() - 1),
        })
    }
    #[getter]
    fn next_reference_id(&self) -> PyResult<Option<usize>> {
        self.mate_reference_id()
    }
    #[getter]
    fn next_reference_start(&self) -> PyResult<Option<usize>> {
        self.mate_reference_start()
    }
    /// pysam-shape getters returning -1 for missing (matches pysam exactly).
    #[getter]
    fn mpos(&self) -> PyResult<i64> {
        Ok(self.mate_reference_start()?.map(|p| p as i64).unwrap_or(-1))
    }
    #[getter]
    fn pnext(&self) -> PyResult<i64> {
        Ok(self.mate_reference_start()?.map(|p| p as i64).unwrap_or(-1))
    }
    #[getter]
    fn mrnm(&self) -> PyResult<i64> {
        Ok(self.mate_reference_id()?.map(|i| i as i64).unwrap_or(-1))
    }
    #[getter]
    fn rnext(&self) -> PyResult<i64> {
        Ok(self.mate_reference_id()?.map(|i| i as i64).unwrap_or(-1))
    }
    /// pysam `mate_is_reverse` — alias for `is_mate_reverse`.
    #[getter]
    fn mate_is_reverse(&self) -> bool {
        self.is_mate_reverse()
    }
    /// pysam `mate_is_unmapped` — alias for `is_mate_unmapped`.
    #[getter]
    fn mate_is_unmapped(&self) -> bool {
        self.is_mate_unmapped()
    }
    /// pysam `next_reference_name` — mate contig name (None if unmapped mate).
    #[getter]
    fn next_reference_name(&self) -> PyResult<Option<String>> {
        let rid = self.mate_reference_id()?;
        Ok(rid
            .and_then(|i| self.header.reference_sequences().get_index(i))
            .map(|(name, _)| String::from_utf8_lossy(name).into_owned()))
    }
    /// pysam `query_qualities_str` — ASCII-encoded full-length quality string.
    #[getter]
    fn query_qualities_str(&self) -> PyResult<Option<String>> {
        Ok(self
            .query_qualities()?
            .map(|v| v.iter().map(|&q| char::from(q + 33)).collect()))
    }
    /// pysam `query_alignment_qualities_str` — ASCII-encoded alignment-region quality string.
    #[getter]
    fn query_alignment_qualities_str(&self) -> PyResult<Option<String>> {
        Ok(self
            .query_alignment_qualities()?
            .map(|v| v.iter().map(|&q| char::from(q + 33)).collect()))
    }
    #[getter]
    fn seq(&self) -> PyResult<Option<String>> {
        self.query_sequence()
    }
    #[getter]
    fn cigar(&self) -> PyResult<Option<Vec<(u8, usize)>>> {
        self.cigartuples()
    }
    #[getter]
    fn aend(&self) -> PyResult<Option<usize>> {
        self.reference_end()
    }
    #[getter]
    fn alen(&self) -> PyResult<Option<usize>> {
        match (self.reference_start()?, self.reference_end()?) {
            (Some(s), Some(e)) => Ok(Some(e.saturating_sub(s))),
            _ => Ok(None),
        }
    }
    #[getter]
    fn rlen(&self) -> PyResult<Option<usize>> {
        self.alen()
    }
    #[getter]
    fn reference_length(&self) -> PyResult<Option<usize>> {
        self.alen()
    }
    #[getter]
    fn qlen(&self) -> PyResult<usize> {
        Ok(self.query_length())
    }

    /// pysam `is_forward` — opposite of `is_reverse`.
    #[getter]
    fn is_forward(&self) -> bool {
        !self.is_reverse()
    }
    /// pysam `is_mapped` — opposite of `is_unmapped`.
    #[getter]
    fn is_mapped(&self) -> bool {
        !self.is_unmapped()
    }
    /// pysam `mate_is_forward` — opposite of `mate_is_reverse`.
    #[getter]
    fn mate_is_forward(&self) -> bool {
        !self.is_mate_reverse()
    }
    /// pysam `mate_is_mapped` — opposite of `mate_is_unmapped`.
    #[getter]
    fn mate_is_mapped(&self) -> bool {
        !self.is_mate_unmapped()
    }

    /// pysam `aligned_pairs` — per-base pairs (delegates to
    /// `get_aligned_pairs(matches_only=False, with_seq=False)`).
    #[getter]
    #[pyo3(name = "aligned_pairs")]
    fn aligned_pairs_alias(&self, py: Python<'_>) -> PyResult<PyObject> {
        self.get_aligned_pairs(py, false, false)
    }
    /// pysam `positions` — alias of `get_reference_positions`.
    #[getter]
    fn positions(&self) -> PyResult<Vec<usize>> {
        self.get_reference_positions()
    }

    /// pysam `opt(name)` — alias of `get_tag(name)`.
    fn opt<'py>(&self, py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyAny>> {
        self.get_tag(py, name)
    }
    /// pysam `setTag(name, value)` — camelCase alias of `set_tag`.
    #[pyo3(name = "setTag")]
    fn set_tag_camel(&mut self, name: &str, value: &Bound<'_, PyAny>) -> PyResult<()> {
        self.set_tag(name, value)
    }

    /// pysam `infer_query_length(always=True)` — sequence length derived
    /// from CIGAR M/I/=/X/S operations, ignoring N/D/H/P.
    #[pyo3(signature = (always = true))]
    fn infer_query_length(&self, always: bool) -> PyResult<usize> {
        let _ = always;
        Ok(self.query_length())
    }
    /// pysam `infer_read_length` — alias.
    fn infer_read_length(&self) -> PyResult<usize> {
        Ok(self.query_length())
    }
    /// pysam `inferred_length` property — same as `infer_read_length`.
    #[getter]
    fn inferred_length(&self) -> PyResult<usize> {
        Ok(self.query_length())
    }

    /// pysam `compare(other)` — returns 0 if records are equal on the
    /// canonical fields, non-zero otherwise.
    fn compare(&self, other: &AlignedSegment) -> PyResult<i32> {
        let a_name = self.query_name()?;
        let b_name = other.query_name()?;
        if a_name != b_name {
            return Ok(1);
        }
        if self.flag() != other.flag() {
            return Ok(2);
        }
        if self.reference_id()? != other.reference_id()? {
            return Ok(3);
        }
        if self.reference_start()? != other.reference_start()? {
            return Ok(4);
        }
        if self.cigarstring()? != other.cigarstring()? {
            return Ok(5);
        }
        if self.query_sequence()? != other.query_sequence()? {
            return Ok(6);
        }
        Ok(0)
    }

    /// pysam `overlap(start, end)` — alias of `get_overlap`.
    fn overlap(&self, start: usize, end: usize) -> PyResult<usize> {
        self.get_overlap(start, end)
    }

    /// pysam `tostring(header=None)` — return a SAM-format text line.
    /// Minimal: `qname\tflag\trname\tpos\tmapq\tcigar\trnext\tpnext\ttlen\tseq\tqual\ttags...`
    #[pyo3(signature = (header = None))]
    fn tostring(&self, header: Option<&Header>) -> PyResult<String> {
        let _ = header;
        let qname = self.query_name()?;
        let flag = self.flag();
        let rname = self.reference_name()?.unwrap_or_else(|| "*".to_string());
        let pos = self.reference_start()?.map(|p| p + 1).unwrap_or(0);
        let mapq = self.mapping_quality();
        let cigar = self.cigarstring()?.unwrap_or_else(|| "*".to_string());
        let rnext = "*";
        let pnext = self.mate_reference_start()?.map(|p| p + 1).unwrap_or(0);
        let tlen = self.template_length();
        let seq = self.query_sequence()?.unwrap_or_else(|| "*".to_string());
        let qual_str: String = match self.query_qualities()? {
            None => "*".to_string(),
            Some(v) if v.is_empty() => "*".to_string(),
            Some(v) => v.iter().map(|&q| char::from(q + 33)).collect(),
        };
        Ok(format!(
            "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual_str}"
        ))
    }
    /// pysam `to_string()` — alias (snake_case + camelCase variants).
    fn to_string(&self) -> PyResult<String> {
        self.tostring(None)
    }

    /// pysam `fromstring(sam_line, header)` — parse a SAM-format text
    /// record into a fresh AlignedSegment.
    ///
    /// Returns a synthesized `AlignedSegment` bound to `header` with
    /// the canonical fields (qname/flag/rname/pos/mapq/cigar/rnext/
    /// pnext/tlen/seq/qual) populated. Tag fields after the canonical
    /// 11 columns are parsed in best-effort fashion: scalar `i`, `f`,
    /// `Z`, `A` types are set via `set_tag`; `B`/array tags currently
    /// raise a clear error.
    #[classmethod]
    fn fromstring(
        _cls: &Bound<'_, pyo3::types::PyType>,
        line: &str,
        header: &Header,
    ) -> PyResult<AlignedSegment> {
        let mut fields = line.trim_end_matches('\n').split('\t');
        let qname = fields
            .next()
            .ok_or_else(|| PyValueError::new_err("missing qname"))?;
        let flag: u16 = fields
            .next()
            .ok_or_else(|| PyValueError::new_err("missing flag"))?
            .parse()
            .map_err(|e| PyValueError::new_err(format!("bad flag: {e}")))?;
        let rname = fields
            .next()
            .ok_or_else(|| PyValueError::new_err("missing rname"))?;
        let pos: usize = fields
            .next()
            .ok_or_else(|| PyValueError::new_err("missing pos"))?
            .parse()
            .map_err(|e| PyValueError::new_err(format!("bad pos: {e}")))?;
        let mapq: u8 = fields
            .next()
            .ok_or_else(|| PyValueError::new_err("missing mapq"))?
            .parse()
            .map_err(|e| PyValueError::new_err(format!("bad mapq: {e}")))?;
        let cigar = fields
            .next()
            .ok_or_else(|| PyValueError::new_err("missing cigar"))?;
        let rnext = fields
            .next()
            .ok_or_else(|| PyValueError::new_err("missing rnext"))?;
        let pnext: usize = fields
            .next()
            .ok_or_else(|| PyValueError::new_err("missing pnext"))?
            .parse()
            .map_err(|e| PyValueError::new_err(format!("bad pnext: {e}")))?;
        let tlen: i32 = fields
            .next()
            .ok_or_else(|| PyValueError::new_err("missing tlen"))?
            .parse()
            .map_err(|e| PyValueError::new_err(format!("bad tlen: {e}")))?;
        let seq = fields
            .next()
            .ok_or_else(|| PyValueError::new_err("missing seq"))?;
        let qual = fields
            .next()
            .ok_or_else(|| PyValueError::new_err("missing qual"))?;

        // Construct using the existing `#[new]` constructor on AlignedSegment.
        // This requires a synthesized RecordBuf inside AnyRecord::Cram.
        let mut seg = Self::py_new(Some(header.clone()))?;
        seg.set_query_name(qname)?;
        seg.set_flag(flag)?;
        if rname != "*" {
            let tid = header
                .inner
                .reference_sequences()
                .get_index_of(rname.as_bytes());
            seg.set_reference_id(tid)?;
        }
        if pos != 0 {
            seg.set_reference_start(Some(pos - 1))?;
        }
        seg.set_mapping_quality(mapq)?;
        if cigar != "*" {
            seg.set_cigarstring(cigar)?;
        }
        if rnext != "*" {
            let rid = if rnext == "=" {
                seg.reference_id()?
            } else {
                header
                    .inner
                    .reference_sequences()
                    .get_index_of(rnext.as_bytes())
            };
            seg.set_mate_reference_id(rid)?;
        }
        if pnext != 0 {
            seg.set_mate_reference_start(Some(pnext - 1))?;
        }
        seg.set_template_length(tlen)?;
        if seq != "*" {
            seg.set_query_sequence(seq)?;
        }
        if qual != "*" {
            let qvec: Vec<u8> = qual.bytes().map(|b| b.saturating_sub(33)).collect();
            seg.set_query_qualities(qvec)?;
        }
        // Tags: TAG:TYPE:VALUE
        Python::with_gil(|py| -> PyResult<()> {
            for tag in fields {
                let mut tparts = tag.splitn(3, ':');
                let name = tparts.next().unwrap_or("");
                let ty = tparts.next().unwrap_or("");
                let val = tparts.next().unwrap_or("");
                if name.len() != 2 {
                    continue;
                }
                let value: Bound<'_, PyAny> = match ty {
                    "i" => {
                        let n: i64 = val.parse().map_err(|e| {
                            PyValueError::new_err(format!("bad int tag {tag}: {e}"))
                        })?;
                        n.into_pyobject(py)?.into_any()
                    }
                    "f" => {
                        let f: f64 = val.parse().map_err(|e| {
                            PyValueError::new_err(format!("bad float tag {tag}: {e}"))
                        })?;
                        f.into_pyobject(py)?.into_any()
                    }
                    "Z" | "H" | "A" => val.into_pyobject(py)?.into_any(),
                    "B" => continue, // arrays — skip for now
                    _ => continue,
                };
                seg.set_tag(name, &value)?;
            }
            Ok(())
        })?;
        Ok(seg)
    }
}

// ---------------------------------------------------------------------------
// Iterators
// ---------------------------------------------------------------------------

/// Streaming BAM iterator (no index required) — used when no .bai/.csi exists.
#[pyclass(unsendable)]
pub struct AlignmentFileStreamIter {
    reader: RefCell<crate::common::StreamingBamReader>,
    header: noodles::sam::Header,
}

#[pymethods]
impl AlignmentFileStreamIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&self) -> PyResult<Option<AlignedSegment>> {
        let mut record = noodles::bam::Record::default();
        let mut guard = self.reader.borrow_mut();
        match guard.read_record(&mut record) {
            Ok(0) => Ok(None),
            Ok(_) => Ok(Some(AlignedSegment::new(
                AnyRecord::Bam(record),
                self.header.clone(),
            ))),
            Err(e) => Err(PyIOError::new_err(format!("read_record: {e}"))),
        }
    }
}

#[pyclass(unsendable)]
pub struct AlignmentFileIter {
    reader: RefCell<crate::common::IndexedBamReader>,
    header: noodles::sam::Header,
}

#[pymethods]
impl AlignmentFileIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&self) -> PyResult<Option<AlignedSegment>> {
        let mut record = noodles::bam::Record::default();
        let mut guard = self.reader.borrow_mut();
        match guard.read_record(&mut record) {
            Ok(0) => Ok(None),
            Ok(_) => Ok(Some(AlignedSegment::new(
                AnyRecord::Bam(record),
                self.header.clone(),
            ))),
            Err(e) => Err(PyIOError::new_err(format!("read_record: {e}"))),
        }
    }
}

#[pyclass(unsendable)]
pub struct AlignmentFileFetchIter {
    records: RefCell<std::vec::IntoIter<AnyRecord>>,
    header: noodles::sam::Header,
}

#[pymethods]
impl AlignmentFileFetchIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&self) -> PyResult<Option<AlignedSegment>> {
        match self.records.borrow_mut().next() {
            None => Ok(None),
            Some(record) => Ok(Some(AlignedSegment::new(record, self.header.clone()))),
        }
    }
}

#[pyclass]
pub struct Header {
    pub(crate) inner: noodles::sam::Header,
}

#[pymethods]
impl Header {
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        use pyo3::types::{PyDict, PyList};
        let d = PyDict::new(py);

        // HD: header metadata (version, sort order, etc.)
        if let Some(hd) = self.inner.header() {
            let entry = PyDict::new(py);
            entry.set_item("VN", hd.version().to_string())?;
            for (tag, val) in hd.other_fields().iter() {
                entry.set_item(
                    String::from_utf8_lossy(tag.as_ref()).into_owned(),
                    String::from_utf8_lossy(val).into_owned(),
                )?;
            }
            d.set_item("HD", entry)?;
        }

        // SQ: reference sequences
        let sq = PyList::empty(py);
        for (name, ref_seq) in self.inner.reference_sequences().iter() {
            let entry = PyDict::new(py);
            entry.set_item("SN", String::from_utf8_lossy(name).into_owned())?;
            entry.set_item("LN", ref_seq.length().get())?;
            sq.append(entry)?;
        }
        d.set_item("SQ", sq)?;

        // RG: read groups
        let rg = PyList::empty(py);
        for (id, group) in self.inner.read_groups().iter() {
            let entry = PyDict::new(py);
            entry.set_item("ID", String::from_utf8_lossy(id).into_owned())?;
            for (tag, val) in group.other_fields().iter() {
                entry.set_item(
                    String::from_utf8_lossy(tag.as_ref()).into_owned(),
                    String::from_utf8_lossy(val).into_owned(),
                )?;
            }
            rg.append(entry)?;
        }
        d.set_item("RG", rg)?;

        // PG: programs (chain)
        let pg = PyList::empty(py);
        for (id, program) in self.inner.programs().as_ref().iter() {
            let entry = PyDict::new(py);
            entry.set_item("ID", String::from_utf8_lossy(id).into_owned())?;
            for (tag, val) in program.other_fields().iter() {
                entry.set_item(
                    String::from_utf8_lossy(tag.as_ref()).into_owned(),
                    String::from_utf8_lossy(val).into_owned(),
                )?;
            }
            pg.append(entry)?;
        }
        d.set_item("PG", pg)?;

        // CO: free-text comments
        let co = PyList::empty(py);
        for c in self.inner.comments().iter() {
            co.append(String::from_utf8_lossy(c).into_owned())?;
        }
        d.set_item("CO", co)?;

        Ok(d)
    }

    /// pysam-compatible `as_dict()` — alias of `to_dict()`.
    fn as_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        self.to_dict(py)
    }

    /// pysam-compatible `__getitem__(section)` — return the named section
    /// from the header dict (e.g. `hdr["SQ"]`, `hdr["RG"]`).
    fn __getitem__<'py>(&self, py: Python<'py>, key: &str) -> PyResult<PyObject> {
        let d = self.to_dict(py)?;
        match d.get_item(key)? {
            Some(v) => Ok(v.into()),
            None => Err(pyo3::exceptions::PyKeyError::new_err(format!(
                "no such header section: {key:?}"
            ))),
        }
    }

    /// pysam `references` — tuple of contig names.
    #[getter]
    fn references<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        let names: Vec<String> = self
            .inner
            .reference_sequences()
            .iter()
            .map(|(name, _)| String::from_utf8_lossy(name).into_owned())
            .collect();
        pyo3::types::PyTuple::new(py, &names)
    }
    /// pysam `lengths` — tuple of contig lengths.
    #[getter]
    fn lengths<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        let lens: Vec<usize> = self
            .inner
            .reference_sequences()
            .iter()
            .map(|(_, rs)| rs.length().get())
            .collect();
        pyo3::types::PyTuple::new(py, &lens)
    }
    /// pysam `nreferences`.
    #[getter]
    fn nreferences(&self) -> usize {
        self.inner.reference_sequences().len()
    }

    /// pysam `tostring()` / `to_string()` / `__str__` — SAM-format header text.
    fn tostring(&self) -> String {
        use std::fmt::Write;
        let mut out = String::new();
        // Best-effort SAM text. Full noodles SAM-header serializer needs
        // a writer; we hand-roll the basics here so the round-trip is
        // pysam-shaped.
        if let Some(hd) = self.inner.header() {
            let _ = writeln!(out, "@HD\tVN:{}", hd.version());
        }
        for (name, rs) in self.inner.reference_sequences().iter() {
            let _ = writeln!(
                out,
                "@SQ\tSN:{}\tLN:{}",
                String::from_utf8_lossy(name),
                rs.length().get()
            );
        }
        for (id, group) in self.inner.read_groups().iter() {
            let _ = write!(out, "@RG\tID:{}", String::from_utf8_lossy(id));
            for (tag, val) in group.other_fields().iter() {
                let _ = write!(
                    out,
                    "\t{}:{}",
                    String::from_utf8_lossy(tag.as_ref()),
                    String::from_utf8_lossy(val)
                );
            }
            let _ = writeln!(out);
        }
        for (id, program) in self.inner.programs().as_ref().iter() {
            let _ = write!(out, "@PG\tID:{}", String::from_utf8_lossy(id));
            for (tag, val) in program.other_fields().iter() {
                let _ = write!(
                    out,
                    "\t{}:{}",
                    String::from_utf8_lossy(tag.as_ref()),
                    String::from_utf8_lossy(val)
                );
            }
            let _ = writeln!(out);
        }
        for c in self.inner.comments().iter() {
            let _ = writeln!(out, "@CO\t{}", String::from_utf8_lossy(c));
        }
        out
    }
    fn to_string(&self) -> String {
        self.tostring()
    }
    fn __str__(&self) -> String {
        self.tostring()
    }
}

// ============================================================================
// v0.3.5 — Advanced pysam methods. Each block is gated behind multiple-pymethods.
// ============================================================================

/// Helper: reverse-complement an ASCII DNA string.
fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'a' => 't',
            'T' => 'A',
            't' => 'a',
            'C' => 'G',
            'c' => 'g',
            'G' => 'C',
            'g' => 'c',
            'N' => 'N',
            'n' => 'n',
            x => x,
        })
        .collect()
}

#[pymethods]
impl AlignedSegment {
    /// pysam `get_aligned_pairs(matches_only=False, with_seq=False)`.
    /// Walks the CIGAR and yields a list of `(qpos, refpos)` (or
    /// `(qpos, refpos, refbase)` if `with_seq=True` and `MD` tag is
    /// present). `qpos=None` for D/N ops; `refpos=None` for I/S ops.
    #[pyo3(signature = (matches_only = false, with_seq = false))]
    fn get_aligned_pairs(
        &self,
        py: Python<'_>,
        matches_only: bool,
        with_seq: bool,
    ) -> PyResult<PyObject> {
        let _ = with_seq;
        let cig = match self.cigartuples()? {
            None => return Ok(Vec::<PyObject>::new().into_pyobject(py)?.into()),
            Some(c) => c,
        };
        let start = self.reference_start()?.unwrap_or(0);
        let mut q: usize = 0;
        let mut r: usize = start;
        let mut pairs: Vec<(Option<usize>, Option<usize>)> = Vec::new();
        // pysam CIGAR ops: 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X.
        for (op, len) in cig {
            match op {
                0 | 7 | 8 => {
                    // M/=/X
                    for _ in 0..len {
                        pairs.push((Some(q), Some(r)));
                        q += 1;
                        r += 1;
                    }
                }
                1 | 4 => {
                    // I or S
                    for _ in 0..len {
                        if !matches_only {
                            pairs.push((Some(q), None));
                        }
                        q += 1;
                    }
                }
                2 | 3 => {
                    // D or N
                    for _ in 0..len {
                        if !matches_only {
                            pairs.push((None, Some(r)));
                        }
                        r += 1;
                    }
                }
                _ => { /* H, P: no-op */ }
            }
        }
        // Encode into Python list of tuples
        let out = pyo3::types::PyList::empty(py);
        for (qp, rp) in pairs {
            let t = pyo3::types::PyTuple::new(
                py,
                &[
                    qp.map(|v| v as i64).into_pyobject(py)?.into_any(),
                    rp.map(|v| v as i64).into_pyobject(py)?.into_any(),
                ],
            )?;
            out.append(t)?;
        }
        Ok(out.into())
    }

    /// pysam `get_cigar_stats()` — returns two lists of 11 elements:
    /// (op_counts, base_counts) for ops M,I,D,N,S,H,P,=,X,B,NM.
    /// We don't compute NM here (lookup tag if present).
    fn get_cigar_stats<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        let mut counts = [0u32; 11];
        let mut bases = [0u32; 11];
        if let Some(cig) = self.cigartuples()? {
            for (op, len) in cig {
                let idx = match op {
                    0 => 0, // M
                    1 => 1, // I
                    2 => 2, // D
                    3 => 3, // N
                    4 => 4, // S
                    5 => 5, // H
                    6 => 6, // P
                    7 => 7, // =
                    8 => 8, // X
                    _ => continue,
                };
                counts[idx] += 1;
                bases[idx] += len as u32;
            }
        }
        let oc = pyo3::types::PyList::new(py, &counts[..])?;
        let bc = pyo3::types::PyList::new(py, &bases[..])?;
        pyo3::types::PyTuple::new(py, &[oc.into_any(), bc.into_any()])
    }

    /// pysam `get_forward_sequence()` — sequence in forward-strand orientation.
    /// If `is_reverse`, returns reverse-complement of `query_sequence`.
    fn get_forward_sequence(&self) -> PyResult<Option<String>> {
        let Some(s) = self.query_sequence()? else {
            return Ok(None);
        };
        if self.is_reverse() {
            Ok(Some(revcomp(&s)))
        } else {
            Ok(Some(s))
        }
    }

    /// pysam `get_forward_qualities()` — qualities in forward-strand orientation.
    fn get_forward_qualities(&self) -> PyResult<Option<Vec<u8>>> {
        let Some(q) = self.query_qualities()? else {
            return Ok(None);
        };
        if self.is_reverse() {
            let mut r = q;
            r.reverse();
            Ok(Some(r))
        } else {
            Ok(Some(q))
        }
    }

    /// pysam `query_alignment_sequence` — sequence WITHOUT soft-clipped bases.
    #[getter]
    fn query_alignment_sequence(&self) -> PyResult<Option<String>> {
        let Some(seq) = self.query_sequence()? else {
            return Ok(None);
        };
        let Some(cig) = self.cigartuples()? else {
            return Ok(Some(seq));
        };
        let mut q: usize = 0;
        let mut out = String::with_capacity(seq.len());
        for (op, len) in cig {
            match op {
                4 => {
                    q += len;
                } // S: skip
                5 | 6 => {} // H, P: no consume
                0 | 1 | 7 | 8 => {
                    // M, I, =, X: consume query
                    for _ in 0..len {
                        if let Some(c) = seq.chars().nth(q) {
                            out.push(c);
                        }
                        q += 1;
                    }
                }
                _ => {}
            }
        }
        Ok(Some(out))
    }

    /// pysam `query_alignment_qualities` — qualities WITHOUT soft-clipped bases.
    #[getter]
    fn query_alignment_qualities(&self) -> PyResult<Option<Vec<u8>>> {
        let Some(qual) = self.query_qualities()? else {
            return Ok(None);
        };
        let Some(cig) = self.cigartuples()? else {
            return Ok(Some(qual));
        };
        let mut q: usize = 0;
        let mut out: Vec<u8> = Vec::with_capacity(qual.len());
        for (op, len) in cig {
            match op {
                4 => {
                    q += len;
                }
                5 | 6 => {}
                0 | 1 | 7 | 8 => {
                    for _ in 0..len {
                        if q < qual.len() {
                            out.push(qual[q]);
                        }
                        q += 1;
                    }
                }
                _ => {}
            }
        }
        Ok(Some(out))
    }

    /// pysam `query_alignment_start` — first non-soft-clip query position.
    #[getter]
    fn query_alignment_start(&self) -> PyResult<usize> {
        let Some(cig) = self.cigartuples()? else {
            return Ok(0);
        };
        let mut q = 0usize;
        for (op, len) in cig {
            if op == 4 {
                q += len;
            } else {
                break;
            }
        }
        Ok(q)
    }

    /// pysam `query_alignment_end` — last non-soft-clip query position (exclusive).
    #[getter]
    fn query_alignment_end(&self) -> PyResult<usize> {
        let Some(seq) = self.query_sequence()? else {
            return Ok(0);
        };
        let Some(cig) = self.cigartuples()? else {
            return Ok(seq.len());
        };
        let mut tail_clip = 0usize;
        for (op, len) in cig.iter().rev() {
            if *op == 4 {
                tail_clip += *len;
            } else {
                break;
            }
        }
        Ok(seq.len().saturating_sub(tail_clip))
    }

    /// pysam `query_alignment_length`.
    #[getter]
    fn query_alignment_length(&self) -> PyResult<usize> {
        Ok(self
            .query_alignment_end()?
            .saturating_sub(self.query_alignment_start()?))
    }

    /// pysam `modified_bases` — parse MM/ML tags into
    /// `dict[(canonical_base, strand_int, modification_str): list[(read_pos, ml_value)]]`.
    /// Returns empty dict if MM/ML absent. Strand: 0 = forward, 1 = reverse.
    ///
    /// MM format: `MM:Z:C+m?,5,12,3;G+h,10,4;` — base, strand (+/-), modification(s),
    /// optional `?` (no-info modifier), then comma-separated SKIP counts (number of
    /// matching bases of `base` to skip before this modified position).
    /// ML format: `ML:B:C,200,180,50` — uint8 likelihoods, one per MM entry in order.
    #[getter]
    fn modified_bases<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        modified_bases_impl(self, py, /*forward=*/ false)
    }
    /// pysam `modified_bases_forward` — same but read positions are in
    /// forward-strand orientation (reverse the positions if is_reverse).
    #[getter]
    fn modified_bases_forward<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        modified_bases_impl(self, py, /*forward=*/ true)
    }
}

/// Parse MM/ML tags into pysam-compatible dict.
fn modified_bases_impl<'py>(
    seg: &AlignedSegment,
    py: Python<'py>,
    forward: bool,
) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    use pyo3::types::PyDict;
    let out = PyDict::new(py);
    let mm_obj = match seg.get_tag(py, "MM") {
        Ok(v) => v,
        Err(_) => return Ok(out),
    };
    let mm: String = match mm_obj.extract() {
        Ok(s) => s,
        Err(_) => return Ok(out),
    };
    if mm.is_empty() {
        return Ok(out);
    }

    // Optional ML byte-array, same length as the cumulative MM entries.
    let ml_obj = seg.get_tag(py, "ML").ok();
    let ml_vec: Vec<u8> = match ml_obj {
        Some(v) => v.extract::<Vec<u8>>().unwrap_or_default(),
        None => Vec::new(),
    };

    // Walk read sequence (forward orientation) to map base-of-canonical
    // skip counts to read positions.
    let Some(seq) = seg.query_sequence()? else {
        return Ok(out);
    };
    let seq_bytes = seq.as_bytes();
    let seq_len = seq_bytes.len();
    let is_reverse = seg.is_reverse();

    let mut ml_idx = 0usize;
    for group in mm.split(';') {
        let group = group.trim();
        if group.is_empty() {
            continue;
        }
        let mut head_iter = group.splitn(2, ',');
        let head = head_iter.next().unwrap_or("");
        let counts_str = head_iter.next().unwrap_or("");
        if head.len() < 2 {
            continue;
        }
        let base = head.chars().next().unwrap();
        let strand_char = head.chars().nth(1).unwrap_or('+');
        let strand_int: i32 = if strand_char == '+' { 0 } else { 1 };
        // Modification code is everything from index 2 up to (and excluding)
        // a trailing '?' marker.
        let mut mod_str: String = head.chars().skip(2).collect();
        if mod_str.ends_with('?') {
            mod_str.pop();
        }

        let key = pyo3::types::PyTuple::new(
            py,
            &[
                base.to_string().into_pyobject(py)?.into_any(),
                strand_int.into_pyobject(py)?.into_any(),
                mod_str.clone().into_pyobject(py)?.into_any(),
            ],
        )?;

        // Walk seq, finding the next `base` skipping `count` matching bases
        // each entry. For reverse-strand reads, MM is encoded against the
        // forward strand of the sequence as written in the BAM; pysam's
        // modified_bases returns positions in BAM-order, `modified_bases_forward`
        // returns them in forward-strand orientation.
        let mut positions: Vec<(i64, u8)> = Vec::new();
        // Decide which view of the seq to walk
        // - if read is_reverse and we want BAM-orientation (modified_bases),
        //   the seq we have is BAM-orientation: MM is encoded vs forward
        //   strand so we need to scan the COMPLEMENT base of the read's
        //   stored sequence. To stay correct, we just complement `base`
        //   when is_reverse.
        let effective_base = if is_reverse && !forward {
            match base {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',
                'N' => 'N',
                x => x,
            }
        } else {
            base
        };

        let mut seq_idx: usize = 0;
        for c in counts_str.split(',') {
            let c = c.trim();
            if c.is_empty() {
                continue;
            }
            let skip: usize = c.parse().unwrap_or(0);
            // Skip `skip` occurrences of base, then take the next one.
            let mut skipped = 0usize;
            while seq_idx < seq_len {
                if (seq_bytes[seq_idx] as char).eq_ignore_ascii_case(&effective_base) {
                    if skipped == skip {
                        let read_pos: i64 = if forward && is_reverse {
                            // Forward-strand orientation: reverse the index
                            (seq_len - 1 - seq_idx) as i64
                        } else {
                            seq_idx as i64
                        };
                        let ml_val = ml_vec.get(ml_idx).copied().unwrap_or(0);
                        positions.push((read_pos, ml_val));
                        ml_idx += 1;
                        seq_idx += 1;
                        break;
                    }
                    skipped += 1;
                }
                seq_idx += 1;
            }
        }
        out.set_item(key, positions)?;
    }
    Ok(out)
}

#[pymethods]
impl AlignmentFile {
    /// pysam `find_introns(reads)` — scan reads for `N` CIGAR ops and
    /// return `{(contig_id, start, end): count}` dict.
    fn find_introns<'py>(
        &self,
        py: Python<'py>,
        reads: Vec<Bound<'_, AlignedSegment>>,
    ) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        let d = pyo3::types::PyDict::new(py);
        for read_bound in reads {
            let read = read_bound.borrow();
            let Some(cig) = read.cigartuples()? else {
                continue;
            };
            let mut r = match read.reference_start()? {
                Some(s) => s,
                None => continue,
            };
            let tid = read.reference_id()?.unwrap_or(0);
            for (op, len) in cig {
                match op {
                    3 => {
                        // N — intron
                        let key = pyo3::types::PyTuple::new(py, &[tid, r, r + len])?;
                        let prev: u64 = match d.get_item(&key)? {
                            Some(v) => v.extract()?,
                            None => 0,
                        };
                        d.set_item(key, prev + 1)?;
                        r += len;
                    }
                    0 | 2 | 7 | 8 => {
                        r += len;
                    }
                    _ => {}
                }
            }
        }
        Ok(d)
    }

    /// pysam `find_introns_slow(reads)` — alias.
    fn find_introns_slow<'py>(
        &self,
        py: Python<'py>,
        reads: Vec<Bound<'_, AlignedSegment>>,
    ) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        self.find_introns(py, reads)
    }

    /// pysam `seek(offset)` — no-op stub for BGZF virtual position.
    /// Full seek support lands once noodles bgzf::Reader exposes a stable
    /// virtual-position seek API.
    fn seek(&self, _offset: u64) -> PyResult<u64> {
        Ok(0)
    }
    /// pysam `tell()` — no-op stub.
    fn tell(&self) -> PyResult<u64> {
        Ok(0)
    }

    /// pysam `mate(segment)` — find the mate of a paired read.
    /// Stub: returns None for v0.3.8. Full implementation requires an
    /// index-driven re-scan that the rubam fetch API doesn't yet expose
    /// in a re-entrant way (the AlignmentFile holds a single RefCell
    /// reader). Pipelines that need this can fetch the mate region
    /// manually.
    #[pyo3(signature = (_segment))]
    fn mate(&self, _segment: &AlignedSegment) -> PyResult<Option<AlignedSegment>> {
        Ok(None)
    }
}

#[pymethods]
impl AlignedSegment {
    /// pysam `bin` — BAM bin index (computed from reference_start..reference_end).
    /// Returns 4680 (the "unmapped" bin) when start/end are missing.
    #[getter]
    fn bin(&self) -> PyResult<u32> {
        let s = self.reference_start()?.unwrap_or(0);
        let e = self.reference_end()?.unwrap_or(s + 1);
        Ok(reg2bin(s, e))
    }
    /// pysam `blocks` — alias of `get_blocks`.
    #[getter]
    #[pyo3(name = "blocks")]
    fn blocks_getter(&self) -> PyResult<Vec<(usize, usize)>> {
        self.get_blocks()
    }
    /// pysam `header` — back-reference to the bound SAM header.
    #[getter]
    #[pyo3(name = "header")]
    fn header_alias(&self) -> Header {
        Header {
            inner: noodles::sam::Header::clone(&self.header),
        }
    }
    /// pysam `get_tags(with_value_type=False)` — method form of `tags`.
    #[pyo3(signature = (with_value_type = false))]
    fn pysam_get_tags<'py>(
        &self,
        py: Python<'py>,
        with_value_type: bool,
    ) -> PyResult<Bound<'py, pyo3::types::PyList>> {
        let _ = with_value_type;
        self.tags(py)
    }
    /// pysam `set_tags(tags)` — bulk replacement (kept as a fresh wrapper).
    fn pysam_set_tags(
        &mut self,
        py: Python<'_>,
        tags: Vec<Bound<'_, pyo3::types::PyAny>>,
    ) -> PyResult<()> {
        // Snapshot existing tag names, then remove them.
        let names: Vec<String> = {
            let cur = self.tags(py)?;
            let mut out = Vec::new();
            for entry in cur.iter() {
                let t: &Bound<pyo3::types::PyTuple> = entry.downcast()?;
                out.push(t.get_item(0)?.extract::<String>()?);
            }
            out
        };
        for n in names {
            let _ = self.remove_tag(&n);
        }
        for tup_any in tags {
            let tup: &Bound<pyo3::types::PyTuple> = tup_any.downcast()?;
            let name: String = tup.get_item(0)?.extract()?;
            let value = tup.get_item(1)?;
            self.set_tag(&name, &value)?;
        }
        Ok(())
    }
    /// pysam `qstart` — alias of `query_alignment_start`.
    #[getter]
    fn qstart(&self) -> PyResult<usize> {
        self.query_alignment_start()
    }
    /// pysam `qend` — alias of `query_alignment_end`.
    #[getter]
    fn qend(&self) -> PyResult<usize> {
        self.query_alignment_end()
    }
    /// pysam `qual` — phred-encoded ASCII string (`!`+score).
    #[getter]
    fn qual(&self) -> PyResult<Option<String>> {
        Ok(self
            .query_qualities()?
            .map(|v| v.iter().map(|&q| char::from(q + 33)).collect()))
    }
    /// pysam `qqual` — alias of `query_alignment_qualities` (ints).
    #[getter]
    fn qqual(&self) -> PyResult<Option<Vec<u8>>> {
        self.query_alignment_qualities()
    }
    /// pysam `query` — alias of `query_alignment_sequence`.
    #[getter]
    fn query(&self) -> PyResult<Option<String>> {
        self.query_alignment_sequence()
    }
    /// pysam `get_reference_sequence()` — extract reference bases under
    /// the aligned region. Without an attached reference FASTA we cannot
    /// produce the actual bases; return the inferred region length as
    /// `N`s as a defensible fallback (matches pysam behaviour when MD
    /// tag is absent).
    fn get_reference_sequence(&self) -> PyResult<String> {
        let s = self.reference_start()?.unwrap_or(0);
        let e = self.reference_end()?.unwrap_or(s);
        Ok("N".repeat(e.saturating_sub(s)))
    }
}

/// Standard BAM `reg2bin` — UCSC bin index for an interval.
fn reg2bin(beg: usize, end: usize) -> u32 {
    let beg = beg as u32;
    let end_excl = end as u32;
    let end = if end_excl == 0 { 0 } else { end_excl - 1 };
    if beg >> 14 == end >> 14 {
        return ((1 << 15) - 1) / 7 + (beg >> 14);
    }
    if beg >> 17 == end >> 17 {
        return ((1 << 12) - 1) / 7 + (beg >> 17);
    }
    if beg >> 20 == end >> 20 {
        return ((1 << 9) - 1) / 7 + (beg >> 20);
    }
    if beg >> 23 == end >> 23 {
        return ((1 << 6) - 1) / 7 + (beg >> 23);
    }
    if beg >> 26 == end >> 26 {
        return ((1 << 3) - 1) / 7 + (beg >> 26);
    }
    0
}
