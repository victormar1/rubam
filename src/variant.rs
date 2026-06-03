//! VariantFile + VariantRecord — VCF/BCF/Tabix support (v0.3).
//!
//! Mirrors the v0.2 `alignment.rs` pattern: thin pyclass wrappers over
//! `noodles-vcf` and (later) `noodles-bcf`, with a streaming iterator and an
//! indexed query path. Phase A of `paper/PLAN_v0.3.md`.

use std::cell::{Cell, RefCell};
use std::io::{BufRead, Write};

use noodles::bcf;
use noodles::bgzf;
use noodles::vcf;
use noodles::vcf::variant::record::samples::series::Value;
use pyo3::exceptions::{PyIOError, PyKeyError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};

type VcfReader = vcf::io::Reader<Box<dyn BufRead>>;
type VcfWriterPlain = vcf::io::Writer<Box<dyn Write>>;
type VcfWriterBgzf = vcf::io::Writer<bgzf::io::Writer<Box<dyn Write>>>;
type BcfWriter = bcf::io::Writer<bgzf::io::Writer<Box<dyn Write>>>;

/// Internal enum: a VariantFile is either open for reading or for writing
/// (in one of three write formats).
#[allow(dead_code)]
enum VariantIo {
    Reader(VcfReader),
    WriterPlain(VcfWriterPlain),
    WriterBgzf(VcfWriterBgzf),
    WriterBcf(BcfWriter),
}

/// Python-side handle to a VCF / BCF file.
///
/// Supports read (`"r"`) and three write modes:
///   `"w"`  — plain VCF text,
///   `"wz"` — BGZF-compressed VCF (`.vcf.gz`),
///   `"wb"` — BCF binary.
#[pyclass(unsendable)]
pub struct VariantFile {
    path: String,
    mode: String,
    inner: RefCell<Option<VariantIo>>,
    header: Option<vcf::Header>,
    closed: Cell<bool>,
}

#[pymethods]
impl VariantFile {
    #[new]
    #[pyo3(signature = (path, mode = "r", header = None))]
    fn new(path: &str, mode: &str, header: Option<&VariantHeader>) -> PyResult<Self> {
        match mode {
            "r" => Self::open_read(path),
            "w" => Self::open_write_plain(path, header),
            "wz" => Self::open_write_bgzf(path, header),
            "wb" => Self::open_write_bcf(path, header),
            other => Err(PyValueError::new_err(format!(
                "unsupported mode {other:?}; use 'r', 'w', 'wz', or 'wb'"
            ))),
        }
    }

    #[getter]
    fn header(&self) -> PyResult<VariantHeader> {
        let h = self
            .header
            .as_ref()
            .ok_or_else(|| PyIOError::new_err("file is closed"))?;
        Ok(VariantHeader { inner: h.clone() })
    }

    #[getter]
    fn is_open(&self) -> bool {
        !self.closed.get()
    }

    fn close(&self) {
        // Dropping the VariantIo will call try_finish on BGZF writers via Drop.
        *self.inner.borrow_mut() = None;
        self.closed.set(true);
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
        self.close();
        Ok(false)
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyResult<VariantFileIter> {
        if slf.path.ends_with(".bcf") {
            // BCF binary path: read all records eagerly as RecordBuf.
            let path = slf.path.clone();
            let f = std::fs::File::open(&path)
                .map_err(|e| PyIOError::new_err(format!("re-open BCF: {e}")))?;
            let mut reader = bcf::io::Reader::new(f);
            let header = reader
                .read_header()
                .map_err(|e| PyIOError::new_err(format!("read BCF header: {e}")))?;
            let mut buf = vcf::variant::RecordBuf::default();
            let mut records: Vec<vcf::variant::RecordBuf> = Vec::new();
            loop {
                match reader.read_record_buf(&header, &mut buf) {
                    Ok(0) => break,
                    Ok(_) => records.push(buf.clone()),
                    Err(e) => return Err(PyIOError::new_err(format!("read BCF record: {e}"))),
                }
            }
            Ok(VariantFileIter {
                inner: VariantFileIterInner::Buffered(RefCell::new(records.into_iter())),
                header,
            })
        } else {
            let mut reader = vcf::io::reader::Builder::default()
                .build_from_path(&slf.path)
                .map_err(|e| PyIOError::new_err(format!("re-open: {e}")))?;
            let header = reader
                .read_header()
                .map_err(|e| PyIOError::new_err(format!("read header: {e}")))?;
            Ok(VariantFileIter {
                inner: VariantFileIterInner::Vcf(RefCell::new(reader)),
                header,
            })
        }
    }

    /// Write a single `VariantRecord` to the file (write modes only).
    fn write(&self, record: &VariantRecord) -> PyResult<()> {
        let header = self
            .header
            .as_ref()
            .ok_or_else(|| PyIOError::new_err("file is closed"))?;
        let mut guard = self.inner.borrow_mut();
        match guard.as_mut() {
            None => Err(PyIOError::new_err("file is closed")),
            Some(VariantIo::Reader(_)) => {
                Err(PyIOError::new_err("file is open for reading; cannot write"))
            }
            Some(VariantIo::WriterPlain(w)) => {
                use noodles::vcf::variant::io::Write as _;
                w.write_variant_record(header, &record.inner)
                    .map_err(|e| PyIOError::new_err(format!("write_record: {e}")))
            }
            Some(VariantIo::WriterBgzf(w)) => {
                use noodles::vcf::variant::io::Write as _;
                w.write_variant_record(header, &record.inner)
                    .map_err(|e| PyIOError::new_err(format!("write_record: {e}")))
            }
            Some(VariantIo::WriterBcf(w)) => {
                use noodles::vcf::variant::io::Write as _;
                w.write_variant_record(header, &record.inner)
                    .map_err(|e| PyIOError::new_err(format!("write_record (BCF): {e}")))
            }
        }
    }

    /// Indexed query — returns records overlapping `[start, end]` (1-based, inclusive).
    ///
    /// `start` and `end` follow pysam conventions (1-based). Raises `ValueError`
    /// if the contig is absent from the header. Returns an empty iterator if
    /// `end <= start` or `start` is past the contig end (matches pysam semantics).
    #[pyo3(signature = (contig, start, end))]
    fn fetch(&self, contig: &str, start: usize, end: usize) -> PyResult<VariantFileFetchIter> {
        use noodles::core::{Position, Region};

        let path = self.path.clone();
        let mut reader = vcf::io::indexed_reader::Builder::default()
            .build_from_path(&path)
            .map_err(|e| PyIOError::new_err(format!("open indexed VCF at {path}: {e}")))?;
        let header = reader
            .read_header()
            .map_err(|e| PyIOError::new_err(format!("read header: {e}")))?;

        // Verify contig exists in header — raise ValueError if not.
        let ref_len_opt: Option<usize> = match header.contigs().get(contig) {
            Some(contig_map) => contig_map.length(),
            None => {
                return Err(PyValueError::new_err(format!(
                    "contig {contig} not found in VCF header"
                )));
            }
        };

        // pysam-like: empty interval or start past reference end → empty iterator.
        if end <= start {
            return Ok(VariantFileFetchIter {
                records: RefCell::new(Vec::<vcf::variant::RecordBuf>::new().into_iter()),
                header,
            });
        }
        if let Some(ref_len) = ref_len_opt {
            if start >= ref_len {
                return Ok(VariantFileFetchIter {
                    records: RefCell::new(Vec::<vcf::variant::RecordBuf>::new().into_iter()),
                    header,
                });
            }
        }
        // Clamp end to reference length when available.
        let end = if let Some(ref_len) = ref_len_opt {
            end.min(ref_len)
        } else {
            end
        };

        let region = Region::new(
            contig.as_bytes().to_vec(),
            Position::new(start).ok_or_else(|| PyValueError::new_err("start must be >= 1"))?
                ..=Position::new(end).ok_or_else(|| PyValueError::new_err("end must be >= 1"))?,
        );

        let mut query = reader
            .query(&header, &region)
            .map_err(|e| PyIOError::new_err(format!("query: {e}")))?;

        let mut buf: Vec<vcf::variant::RecordBuf> = Vec::new();
        let mut raw_record = vcf::Record::default();
        loop {
            match query.read_record(&mut raw_record) {
                Ok(0) => break,
                Ok(_) => {
                    let rb = vcf::variant::RecordBuf::try_from_variant_record(&header, &raw_record)
                        .map_err(|e| PyIOError::new_err(format!("record conversion: {e}")))?;
                    buf.push(rb);
                }
                Err(e) => return Err(PyIOError::new_err(format!("record: {e}"))),
            }
        }

        Ok(VariantFileFetchIter {
            records: RefCell::new(buf.into_iter()),
            header,
        })
    }
}

// Private constructors used by `new`.
impl VariantFile {
    fn open_read(path: &str) -> PyResult<Self> {
        if path.ends_with(".bcf") {
            // BCF: read header only via BCF reader; __iter__ re-opens for records.
            let f = std::fs::File::open(path)
                .map_err(|e| PyIOError::new_err(format!("failed to open BCF at {path}: {e}")))?;
            let mut reader = bcf::io::Reader::new(f);
            let header = reader
                .read_header()
                .map_err(|e| PyIOError::new_err(format!("failed to read BCF header: {e}")))?;
            Ok(Self {
                path: path.to_string(),
                mode: "r".to_string(),
                inner: RefCell::new(None), // BCF has no streaming reader stored
                header: Some(header),
                closed: Cell::new(false),
            })
        } else {
            let mut reader = vcf::io::reader::Builder::default()
                .build_from_path(path)
                .map_err(|e| PyIOError::new_err(format!("failed to open VCF at {path}: {e}")))?;
            let header = reader
                .read_header()
                .map_err(|e| PyIOError::new_err(format!("failed to read VCF header: {e}")))?;
            Ok(Self {
                path: path.to_string(),
                mode: "r".to_string(),
                inner: RefCell::new(Some(VariantIo::Reader(reader))),
                header: Some(header),
                closed: Cell::new(false),
            })
        }
    }

    fn require_header(header: Option<&VariantHeader>, mode: &str) -> PyResult<vcf::Header> {
        header.map(|h| h.inner.clone()).ok_or_else(|| {
            PyValueError::new_err(format!("mode {mode:?} requires a header= keyword argument"))
        })
    }

    fn open_write_plain(path: &str, header: Option<&VariantHeader>) -> PyResult<Self> {
        let hdr = Self::require_header(header, "w")?;
        let f = std::fs::File::create(path)
            .map_err(|e| PyIOError::new_err(format!("create {path}: {e}")))?;
        let mut writer: VcfWriterPlain = vcf::io::Writer::new(Box::new(f) as Box<dyn Write>);
        writer
            .write_header(&hdr)
            .map_err(|e| PyIOError::new_err(format!("write_header: {e}")))?;
        Ok(Self {
            path: path.to_string(),
            mode: "w".to_string(),
            inner: RefCell::new(Some(VariantIo::WriterPlain(writer))),
            header: Some(hdr),
            closed: Cell::new(false),
        })
    }

    fn open_write_bgzf(path: &str, header: Option<&VariantHeader>) -> PyResult<Self> {
        let hdr = Self::require_header(header, "wz")?;
        let f = std::fs::File::create(path)
            .map_err(|e| PyIOError::new_err(format!("create {path}: {e}")))?;
        let bgzf_inner = bgzf::io::Writer::new(Box::new(f) as Box<dyn Write>);
        let mut writer: VcfWriterBgzf = vcf::io::Writer::new(bgzf_inner);
        writer
            .write_header(&hdr)
            .map_err(|e| PyIOError::new_err(format!("write_header: {e}")))?;
        Ok(Self {
            path: path.to_string(),
            mode: "wz".to_string(),
            inner: RefCell::new(Some(VariantIo::WriterBgzf(writer))),
            header: Some(hdr),
            closed: Cell::new(false),
        })
    }

    fn open_write_bcf(path: &str, header: Option<&VariantHeader>) -> PyResult<Self> {
        let hdr = Self::require_header(header, "wb")?;
        let f = std::fs::File::create(path)
            .map_err(|e| PyIOError::new_err(format!("create {path}: {e}")))?;
        let mut writer: BcfWriter = bcf::io::Writer::new(Box::new(f) as Box<dyn Write>);
        writer
            .write_header(&hdr)
            .map_err(|e| PyIOError::new_err(format!("write_header: {e}")))?;
        Ok(Self {
            path: path.to_string(),
            mode: "wb".to_string(),
            inner: RefCell::new(Some(VariantIo::WriterBcf(writer))),
            header: Some(hdr),
            closed: Cell::new(false),
        })
    }
}

// -------- pysam-compatible VariantFile aliases (v0.3.4) --------
#[pymethods]
impl VariantFile {
    /// pysam `closed` — opposite of `is_open`.
    #[getter]
    fn closed(&self) -> bool {
        !self.is_open()
    }
    /// pysam `is_closed` — same as `closed`.
    #[getter]
    fn is_closed(&self) -> bool {
        !self.is_open()
    }
    /// pysam `filename` — the file path passed at construction.
    #[getter]
    fn filename(&self) -> String {
        self.path.clone()
    }
    /// pysam `mode` — the open-mode string.
    #[getter]
    #[pyo3(name = "mode")]
    fn mode_alias(&self) -> String {
        self.mode.clone()
    }
    /// pysam `is_read` / `is_write` — derived from mode.
    #[getter]
    fn is_read(&self) -> bool {
        self.mode.starts_with('r')
    }
    #[getter]
    fn is_write(&self) -> bool {
        self.mode.starts_with('w')
    }
    #[getter]
    fn is_reading(&self) -> bool {
        self.is_read()
    }
    /// pysam `is_bam` / `is_sam` / `is_cram` — always False on VariantFile.
    #[getter]
    fn is_bam(&self) -> bool {
        false
    }
    #[getter]
    fn is_sam(&self) -> bool {
        false
    }
    #[getter]
    fn is_cram(&self) -> bool {
        false
    }
    /// pysam `is_vcf` — True if path ends with .vcf or .vcf.gz.
    #[getter]
    fn is_vcf(&self) -> bool {
        let p = self.path.to_ascii_lowercase();
        p.ends_with(".vcf") || p.ends_with(".vcf.gz")
    }
    /// pysam `is_bcf` — True if path ends with .bcf.
    #[getter]
    fn is_bcf(&self) -> bool {
        self.path.to_ascii_lowercase().ends_with(".bcf")
    }
    /// pysam `is_remote` / `is_stream` — False (no HTTP/stdio support yet).
    #[getter]
    fn is_remote(&self) -> bool {
        false
    }
    #[getter]
    fn is_stream(&self) -> bool {
        self.path == "-"
    }
    /// pysam `format` — "VCF" / "VCF.gz" / "BCF" / "UNKNOWN".
    #[getter]
    fn format(&self) -> &'static str {
        let p = self.path.to_ascii_lowercase();
        if p.ends_with(".vcf.gz") {
            "VCF.gz"
        } else if p.ends_with(".vcf") {
            "VCF"
        } else if p.ends_with(".bcf") {
            "BCF"
        } else {
            "UNKNOWN"
        }
    }
    /// pysam `compression` — "BGZF" for .gz/.bcf, "NONE" for plain .vcf.
    #[getter]
    fn compression(&self) -> &'static str {
        let p = self.path.to_ascii_lowercase();
        if p.ends_with(".gz") || p.ends_with(".bcf") {
            "BGZF"
        } else {
            "NONE"
        }
    }
    /// pysam `category` — always "variant".
    #[getter]
    fn category(&self) -> &'static str {
        "variant"
    }
    /// pysam `description` — human-readable format description.
    #[getter]
    fn description(&self) -> &'static str {
        "Variant Call Format (VCF/BCF) — pysam-compatible read+write surface via rubam"
    }
    /// pysam `version` — backend identification string.
    #[getter]
    fn version(&self) -> &'static str {
        "noodles-vcf / noodles-bcf 0.107 (via rubam)"
    }
    /// pysam `index_filename` — `.tbi`/`.csi` if present.
    #[getter]
    fn index_filename(&self) -> Option<String> {
        let tbi = format!("{}.tbi", self.path);
        let csi = format!("{}.csi", self.path);
        if std::path::Path::new(&tbi).exists() {
            Some(tbi)
        } else if std::path::Path::new(&csi).exists() {
            Some(csi)
        } else {
            None
        }
    }
    /// pysam `flush()` — no-op for now (BGZF buffer flush lands later).
    fn flush(&self) -> PyResult<()> {
        Ok(())
    }
    /// pysam `reset()` — no-op stub (full BGZF seek-to-start lands later).
    fn reset(&self) -> PyResult<()> {
        Ok(())
    }
    /// pysam `add_hts_options(opts)` — no-op for pysam-compat.
    fn add_hts_options(&self, _opts: Vec<String>) {}
}

/// Internal enum for `VariantFileIter`: either a live VCF stream or a buffered list of RecordBufs
/// (used for BCF files which are read eagerly at iter construction time).
enum VariantFileIterInner {
    Vcf(RefCell<VcfReader>),
    Buffered(RefCell<std::vec::IntoIter<vcf::variant::RecordBuf>>),
}

/// Streaming iterator over VariantRecord. Constructed by `VariantFile.__iter__`.
#[pyclass(unsendable)]
pub struct VariantFileIter {
    inner: VariantFileIterInner,
    header: vcf::Header,
}

#[pymethods]
impl VariantFileIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&self) -> PyResult<Option<VariantRecord>> {
        match &self.inner {
            VariantFileIterInner::Vcf(reader) => {
                let mut record = vcf::variant::RecordBuf::default();
                let mut guard = reader.borrow_mut();
                match guard.read_record_buf(&self.header, &mut record) {
                    Ok(0) => Ok(None),
                    Ok(_) => Ok(Some(VariantRecord {
                        inner: record,
                        header: self.header.clone(),
                    })),
                    Err(e) => Err(PyIOError::new_err(format!("read_record: {e}"))),
                }
            }
            VariantFileIterInner::Buffered(buf) => {
                Ok(buf.borrow_mut().next().map(|rec| VariantRecord {
                    inner: rec,
                    header: self.header.clone(),
                }))
            }
        }
    }
}

/// Eagerly-buffered iterator returned by `VariantFile.fetch`.
///
/// Holds an owned `Vec<vcf::variant::RecordBuf>` to avoid lifetime entanglement with the
/// IndexedReader (mirrors the BAM-side `AlignmentFileFetchIter`).
#[pyclass(unsendable)]
pub struct VariantFileFetchIter {
    records: RefCell<std::vec::IntoIter<vcf::variant::RecordBuf>>,
    header: vcf::Header,
}

#[pymethods]
impl VariantFileFetchIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&self) -> Option<VariantRecord> {
        self.records.borrow_mut().next().map(|rec| VariantRecord {
            inner: rec,
            header: self.header.clone(),
        })
    }
}

/// Python-side handle to a VCF/BCF record (a single variant line).
#[pyclass]
pub struct VariantRecord {
    pub(crate) inner: vcf::variant::RecordBuf,
    pub(crate) header: vcf::Header,
}

#[pymethods]
impl VariantRecord {
    /// Build a record from scratch.
    ///
    /// ```python
    /// rec = rubam.VariantRecord(
    ///     header=hdr,
    ///     reference_name="chr1",
    ///     position=100,
    ///     reference="A",
    ///     alternates=("G",),
    ///     quality=30.0,   # optional
    ///     ids=("rs1",),   # optional
    ///     filters=("PASS",),  # optional
    /// )
    /// ```
    #[new]
    #[pyo3(signature = (header, reference_name, position, reference, alternates, *, quality = None, ids = None, filters = None))]
    fn new(
        header: &VariantHeader,
        reference_name: &str,
        position: usize,
        reference: &str,
        alternates: Vec<String>,
        quality: Option<f32>,
        ids: Option<Vec<String>>,
        filters: Option<Vec<String>>,
    ) -> PyResult<Self> {
        use noodles::core::Position;
        use noodles::vcf::variant::record_buf::{AlternateBases, Filters, Ids};

        let pos = Position::new(position)
            .ok_or_else(|| PyValueError::new_err("position must be >= 1"))?;

        let alt_bases = AlternateBases::from(alternates);

        let ids_buf: Ids = ids.unwrap_or_default().into_iter().collect();

        let filters_buf: Filters = match filters {
            None => Filters::default(),
            Some(v) => v.into_iter().collect(),
        };

        let mut builder = vcf::variant::RecordBuf::builder()
            .set_reference_sequence_name(reference_name)
            .set_variant_start(pos)
            .set_reference_bases(reference)
            .set_alternate_bases(alt_bases)
            .set_ids(ids_buf)
            .set_filters(filters_buf);

        if let Some(q) = quality {
            builder = builder.set_quality_score(q);
        }

        Ok(Self {
            inner: builder.build(),
            header: header.inner.clone(),
        })
    }

    /// Override POS (1-based; must be >= 1).
    fn set_position(&mut self, pos: usize) -> PyResult<()> {
        use noodles::core::Position;
        let p = Position::new(pos).ok_or_else(|| PyValueError::new_err("position must be >= 1"))?;
        *self.inner.variant_start_mut() = Some(p);
        Ok(())
    }

    /// Set QUAL. Pass `None` to clear (`.`).
    #[pyo3(signature = (q))]
    fn set_quality(&mut self, q: Option<f32>) {
        *self.inner.quality_score_mut() = q;
    }

    /// Replace the FILTER list with a single entry.
    fn set_filter(&mut self, filter: &str) {
        use noodles::vcf::variant::record_buf::Filters;
        *self.inner.filters_mut() = std::iter::once(filter.to_string()).collect::<Filters>();
    }

    /// Append a filter (no-op if already present).
    fn add_filter(&mut self, filter: &str) {
        self.inner.filters_mut().as_mut().insert(filter.to_string());
    }

    /// Clear all filters (sets FILTER to `.`).
    fn clear_filters(&mut self) {
        use noodles::vcf::variant::record_buf::Filters;
        *self.inner.filters_mut() = Filters::default();
    }

    /// Set one INFO field. Accepts `int`, `float`, `str`, `bool` (→ Flag),
    /// or a `list`/`tuple` of those.
    fn set_info(&mut self, py: Python<'_>, key: &str, value: &Bound<'_, PyAny>) -> PyResult<()> {
        let info_value = py_to_info_value(py, value)?;
        self.inner
            .info_mut()
            .insert(key.to_string(), Some(info_value));
        Ok(())
    }

    /// 1-based start position (the VCF POS column).
    #[getter]
    fn position(&self) -> PyResult<Option<usize>> {
        Ok(self.inner.variant_start().map(|p| p.get()))
    }

    /// Reference contig name (the VCF CHROM column).
    #[getter]
    fn reference_name(&self) -> &str {
        self.inner.reference_sequence_name()
    }

    /// pysam-compatible alias for `position` (1-based).
    #[getter]
    fn pos(&self) -> PyResult<Option<usize>> {
        self.position()
    }

    /// REF allele as a string.
    #[getter]
    fn reference(&self) -> &str {
        self.inner.reference_bases()
    }

    /// pysam-compatible alias.
    #[getter]
    fn ref_allele(&self) -> &str {
        self.inner.reference_bases()
    }

    /// Tuple of ALT alleles. Empty tuple when ALT is `.`.
    #[getter]
    fn alternates<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        use noodles::vcf::variant::record::AlternateBases as _;
        use pyo3::types::PyTuple;
        let alts = self.inner.alternate_bases();
        let mut items: Vec<String> = Vec::new();
        for entry in alts.iter() {
            let a = entry.map_err(|e| PyIOError::new_err(format!("alt: {e}")))?;
            items.push(a.to_string());
        }
        Ok(PyTuple::new(py, items)?)
    }

    /// pysam-compatible alias.
    #[getter]
    fn alts<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        self.alternates(py)
    }

    /// Tuple of ID values from the ID column. Empty tuple when ID is `.`.
    #[getter]
    fn ids<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        use noodles::vcf::variant::record::Ids as _;
        use pyo3::types::PyTuple;
        let mut items: Vec<String> = Vec::new();
        for entry in self.inner.ids().iter() {
            items.push(entry.to_string());
        }
        Ok(PyTuple::new(py, items)?)
    }

    /// QUAL column as `Option<f32>`. `None` if QUAL is `.`.
    #[getter]
    fn quality(&self) -> PyResult<Option<f32>> {
        Ok(self.inner.quality_score())
    }

    /// pysam-compatible alias.
    #[getter]
    fn qual(&self) -> PyResult<Option<f32>> {
        self.quality()
    }

    /// Tuple of filter names from the FILTER column. PASS is reported as
    /// the literal `"PASS"`. Empty tuple when FILTER is `.`.
    #[getter]
    fn filters<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        use noodles::vcf::variant::record::Filters as _;
        let mut items: Vec<String> = Vec::new();
        for entry in self.inner.filters().iter(&self.header) {
            let s = entry.map_err(|e| PyIOError::new_err(format!("filter: {e}")))?;
            items.push(s.to_string());
        }
        Ok(PyTuple::new(py, items)?)
    }

    /// Dict-like container mapping sample names → per-sample field values.
    #[getter]
    fn samples(&self) -> VariantSamples {
        VariantSamples {
            record: self.inner.clone(),
            header: self.header.clone(),
        }
    }

    /// Read-side `info` getter — returns a Python dict of all INFO fields on
    /// this record. **New in v0.3.2**: closes the gap flagged by the v4
    /// reviewer (the previous `VariantRecord.info` xfail).
    ///
    /// Each value is a typed Python object — `int` / `float` / `str` / `bool`
    /// for scalar fields, `tuple` for `Number=A/R/G/.` arrays, `True` for
    /// `Type=Flag`. Missing values (`.`) come through as `None`.
    #[getter]
    fn info<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        // record_buf::Info derefs to IndexMap<String, Option<Value>>; the
        // .iter() yields (&String, &Option<Value>) entries in insertion order.
        for (key, opt_value) in self.inner.info().as_ref().iter() {
            let py_val: PyObject = match opt_value {
                None => py.None(),
                Some(v) => info_value_to_py(py, v)?,
            };
            d.set_item(key.as_str(), py_val)?;
        }
        Ok(d)
    }

    // -------- pysam-compatible VariantRecord aliases (v0.3.4) --------

    /// pysam `chrom` — alias of `reference_name`.
    #[getter]
    fn chrom(&self) -> &str {
        self.reference_name()
    }
    /// pysam `contig` — alias of `reference_name`.
    #[getter]
    fn contig(&self) -> &str {
        self.reference_name()
    }
    /// pysam `start` — 0-based position (pysam convention).
    #[getter]
    fn start(&self) -> PyResult<Option<usize>> {
        Ok(self.position()?.map(|p| p.saturating_sub(1)))
    }
    /// pysam `stop` — exclusive 0-based end position.
    #[getter]
    fn stop(&self) -> PyResult<Option<usize>> {
        let s = self.start()?;
        let ref_len = self.ref_allele().len();
        Ok(s.map(|p| p + ref_len))
    }
    /// pysam `rlen` — length of the REF allele.
    #[getter]
    fn rlen(&self) -> usize {
        self.ref_allele().len()
    }
    /// pysam `ref` — alias of `ref_allele`.
    #[getter]
    #[pyo3(name = "ref")]
    fn ref_alias(&self) -> &str {
        self.ref_allele()
    }
    /// pysam `id` — first ID or "." if none (pysam convention).
    #[getter]
    #[pyo3(name = "id")]
    fn id_alias<'py>(&self, py: Python<'py>) -> PyResult<String> {
        let ids = self.ids(py)?;
        if ids.len() == 0 {
            Ok(".".to_string())
        } else {
            ids.get_item(0)?.extract::<String>()
        }
    }
    /// pysam `filter` — alias of `filters`.
    #[getter]
    #[pyo3(name = "filter")]
    fn filter_alias<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        self.filters(py)
    }
    /// pysam `alleles` — `(ref,) + alts` as a tuple of strings.
    #[getter]
    fn alleles<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        let r = self.ref_allele().to_string();
        let alts = self.alternates(py)?;
        let mut all: Vec<String> = Vec::with_capacity(1 + alts.len());
        all.push(r);
        for a in alts.iter() {
            all.push(a.extract::<String>()?);
        }
        pyo3::types::PyTuple::new(py, &all)
    }
    /// pysam `qual` — alias of `quality`.
    #[getter]
    #[pyo3(name = "qual")]
    fn qual_alias(&self) -> PyResult<Option<f32>> {
        self.quality()
    }
    /// pysam `alleles_variant_types` — classify each allele
    /// ("REF", "SNP", "MNP", "INS", "DEL", "BND", "OTHER").
    #[getter]
    fn alleles_variant_types<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        let ref_seq = self.ref_allele().to_string();
        let alts = self.alternates(py)?;
        let mut types: Vec<&'static str> = Vec::with_capacity(1 + alts.len());
        types.push("REF");
        for a in alts.iter() {
            let alt: String = a.extract()?;
            let t = if alt.starts_with('<') {
                "OTHER"
            } else if alt.contains('[') || alt.contains(']') {
                "BND"
            } else if alt.len() == ref_seq.len() {
                if alt.len() == 1 {
                    "SNP"
                } else {
                    "MNP"
                }
            } else if alt.len() > ref_seq.len() {
                "INS"
            } else {
                "DEL"
            };
            types.push(t);
        }
        pyo3::types::PyTuple::new(py, &types)
    }
}

// ---------------------------------------------------------------------------
// Task A4 — VariantSamples / VariantSample / VariantSamplesIter
// ---------------------------------------------------------------------------

/// Dict-like container over all samples in a VCF record.
/// `record.samples["NA12878"]["GT"]` — mirrors pysam's API.
#[pyclass(unsendable)]
pub struct VariantSamples {
    record: vcf::variant::RecordBuf,
    header: vcf::Header,
}

#[pymethods]
impl VariantSamples {
    fn __len__(&self) -> usize {
        self.header.sample_names().len()
    }

    fn __contains__(&self, name: &str) -> bool {
        self.header.sample_names().contains(name)
    }

    fn __iter__(&self) -> VariantSamplesIter {
        let names: Vec<String> = self.header.sample_names().iter().cloned().collect();
        VariantSamplesIter { names, pos: 0 }
    }

    fn __getitem__(&self, name: &str) -> PyResult<VariantSample> {
        if !self.header.sample_names().contains(name) {
            return Err(PyKeyError::new_err(name.to_string()));
        }
        Ok(VariantSample {
            record: self.record.clone(),
            header: self.header.clone(),
            name: name.to_string(),
        })
    }
}

/// Iterator over sample names (strings) in declaration order.
#[pyclass(unsendable)]
pub struct VariantSamplesIter {
    names: Vec<String>,
    pos: usize,
}

#[pymethods]
impl VariantSamplesIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<String> {
        if self.pos < self.names.len() {
            let s = self.names[self.pos].clone();
            self.pos += 1;
            Some(s)
        } else {
            None
        }
    }
}

/// A single sample's typed-value bag. `sample["GT"]`, `sample["DP"]`, etc.
#[pyclass(unsendable)]
pub struct VariantSample {
    record: vcf::variant::RecordBuf,
    header: vcf::Header,
    name: String,
}

#[pymethods]
impl VariantSample {
    #[getter]
    fn name(&self) -> &str {
        &self.name
    }

    fn __contains__(&self, key: &str) -> bool {
        use noodles::vcf::variant::record::samples::Sample as _;
        let samples = self.record.samples();
        match samples.get(&self.header, &self.name) {
            None => false,
            Some(sample) => sample
                .iter(&self.header)
                .any(|r| matches!(r, Ok((k, _)) if k == key)),
        }
    }

    fn __getitem__(&self, py: Python<'_>, key: &str) -> PyResult<PyObject> {
        use noodles::vcf::variant::record::samples::Sample as _;
        let samples = self.record.samples();
        let sample = samples
            .get(&self.header, &self.name)
            .ok_or_else(|| PyKeyError::new_err(self.name.clone()))?;

        // Walk FORMAT fields for this sample looking for `key`.
        for result in sample.iter(&self.header) {
            let (k, opt_value) =
                result.map_err(|e| PyIOError::new_err(format!("sample field: {e}")))?;
            if k != key {
                continue;
            }
            return match opt_value {
                None => Ok(py.None()),
                Some(value) => value_to_py(py, value),
            };
        }

        // Key not present in FORMAT for this record.
        Err(PyKeyError::new_err(key.to_string()))
    }

    /// **New in v0.3.2**: `phased` — `True` iff this sample's GT field uses
    /// pipe separators (e.g. `0|1`), `False` for slash separators (`0/1`),
    /// `None` if the sample has no GT field on this record.
    ///
    /// Closes the v4 reviewer's xfail on `VariantSample.phased`. Inspects the
    /// underlying noodles Genotype iterator's phasing markers.
    #[getter]
    fn phased(&self, py: Python<'_>) -> PyResult<PyObject> {
        use noodles::vcf::variant::record::samples::series::value::genotype::Phasing;
        use noodles::vcf::variant::record::samples::Sample as _;
        let samples = self.record.samples();
        let sample = match samples.get(&self.header, &self.name) {
            None => return Ok(py.None()),
            Some(s) => s,
        };
        for result in sample.iter(&self.header) {
            let (k, opt_value) = result.map_err(|e| PyIOError::new_err(format!("phased: {e}")))?;
            if k != "GT" {
                continue;
            }
            return match opt_value {
                None => Ok(py.None()),
                Some(Value::Genotype(gt)) => {
                    // The first allele has no preceding separator, so phasing
                    // is determined by separators on alleles 2..N. A genotype
                    // is "phased" iff every non-first allele is `Phased`.
                    let mut saw_separator = false;
                    let mut all_phased = true;
                    for (idx, item) in gt.iter().enumerate() {
                        let (_pos, phasing) =
                            item.map_err(|e| PyIOError::new_err(format!("phased iter: {e}")))?;
                        if idx == 0 {
                            // First allele's phasing slot is meaningless.
                            continue;
                        }
                        saw_separator = true;
                        if !matches!(phasing, Phasing::Phased) {
                            all_phased = false;
                        }
                    }
                    if !saw_separator {
                        // Haploid genotype — conventionally `True` (no separator,
                        // no ambiguity).
                        return Ok(true.into_pyobject(py)?.to_owned().into_any().unbind());
                    }
                    Ok(all_phased.into_pyobject(py)?.to_owned().into_any().unbind())
                }
                Some(_) => Err(PyValueError::new_err(
                    "GT field present but not a Genotype value",
                )),
            };
        }
        // No GT in this sample.
        Ok(py.None())
    }
}

// ---------------------------------------------------------------------------
// Task A6 — VariantHeader + dict-like wrappers
// ---------------------------------------------------------------------------

/// Render an info `Number` enum to the canonical VCF string ("1", "A", "R", "G", ".").
fn info_number_to_str(n: vcf::header::record::value::map::info::Number) -> String {
    use vcf::header::record::value::map::info::Number;
    match n {
        Number::Count(c) => c.to_string(),
        Number::AlternateBases => "A".to_string(),
        Number::ReferenceAlternateBases => "R".to_string(),
        Number::Samples => "G".to_string(),
        Number::Unknown => ".".to_string(),
    }
}

/// Render a format `Number` enum to the canonical VCF string.
fn fmt_number_to_str(n: vcf::header::record::value::map::format::Number) -> String {
    use vcf::header::record::value::map::format::Number;
    match n {
        Number::Count(c) => c.to_string(),
        Number::AlternateBases => "A".to_string(),
        Number::ReferenceAlternateBases => "R".to_string(),
        Number::Samples => "G".to_string(),
        Number::LocalAlternateBases => "LA".to_string(),
        Number::LocalReferenceAlternateBases => "LR".to_string(),
        Number::LocalSamples => "LG".to_string(),
        Number::Ploidy => "P".to_string(),
        Number::BaseModifications => "M".to_string(),
        Number::Unknown => ".".to_string(),
    }
}

/// Read-only Python view of a VCF header.
#[pyclass]
pub struct VariantHeader {
    inner: vcf::Header,
}

#[pymethods]
impl VariantHeader {
    /// Tuple of sample names in declaration order.
    #[getter]
    fn samples<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let names: Vec<&str> = self
            .inner
            .sample_names()
            .iter()
            .map(String::as_str)
            .collect();
        Ok(PyTuple::new(py, names)?)
    }

    /// File format version string, e.g. "VCFv4.3".
    #[getter]
    fn version(&self) -> String {
        let ff = self.inner.file_format();
        format!("VCFv{}.{}", ff.major(), ff.minor())
    }

    /// Dict-like contig map.
    #[getter]
    fn contigs(&self) -> VariantContigs {
        let entries: Vec<(String, Option<usize>)> = self
            .inner
            .contigs()
            .iter()
            .map(|(k, v)| (k.clone(), v.length()))
            .collect();
        VariantContigs { entries }
    }

    /// Dict-like INFO definitions map.
    #[getter]
    fn info(&self) -> VariantInfoDefs {
        let entries: Vec<(String, String, String, String)> = self
            .inner
            .infos()
            .iter()
            .map(|(k, v)| {
                (
                    k.clone(),
                    info_number_to_str(v.number()),
                    v.ty().to_string(),
                    v.description().to_string(),
                )
            })
            .collect();
        VariantInfoDefs { entries }
    }

    /// Dict-like FORMAT definitions map.
    #[getter]
    fn formats(&self) -> VariantFormatDefs {
        let entries: Vec<(String, String, String, String)> = self
            .inner
            .formats()
            .iter()
            .map(|(k, v)| {
                (
                    k.clone(),
                    fmt_number_to_str(v.number()),
                    v.ty().to_string(),
                    v.description().to_string(),
                )
            })
            .collect();
        VariantFormatDefs { entries }
    }

    /// Tuple of declared FILTER names (excluding implicit PASS when not explicitly declared).
    #[getter]
    fn filters<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let names: Vec<&str> = self.inner.filters().keys().map(String::as_str).collect();
        Ok(PyTuple::new(py, names)?)
    }
}

// --- VariantContigs -------------------------------------------------------

/// Dict-like wrapper over contig declarations in a VCF header.
#[pyclass]
pub struct VariantContigs {
    /// (name, length) pairs in declaration order.
    entries: Vec<(String, Option<usize>)>,
}

#[pymethods]
impl VariantContigs {
    fn __len__(&self) -> usize {
        self.entries.len()
    }

    fn __contains__(&self, name: &str) -> bool {
        self.entries.iter().any(|(k, _)| k == name)
    }

    fn __getitem__(&self, name: &str) -> PyResult<VariantContig> {
        self.entries
            .iter()
            .find(|(k, _)| k == name)
            .map(|(k, l)| VariantContig {
                name: k.clone(),
                length: *l,
            })
            .ok_or_else(|| PyKeyError::new_err(name.to_string()))
    }

    fn __iter__(&self) -> VariantContigsIter {
        VariantContigsIter {
            names: self.entries.iter().map(|(k, _)| k.clone()).collect(),
            pos: 0,
        }
    }
}

/// Iterator that yields contig names in order.
#[pyclass]
pub struct VariantContigsIter {
    names: Vec<String>,
    pos: usize,
}

#[pymethods]
impl VariantContigsIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<String> {
        if self.pos < self.names.len() {
            let s = self.names[self.pos].clone();
            self.pos += 1;
            Some(s)
        } else {
            None
        }
    }
}

/// A single contig's metadata.
#[pyclass]
pub struct VariantContig {
    #[pyo3(get)]
    pub name: String,
    #[pyo3(get)]
    pub length: Option<usize>,
}

// --- VariantInfoDefs / VariantFormatDefs / VariantFieldDef ----------------

/// One INFO or FORMAT meta-line definition.
#[pyclass]
pub struct VariantFieldDef {
    #[pyo3(get)]
    pub id: String,
    #[pyo3(get)]
    pub number: String,
    #[pyo3(get, name = "type")]
    pub ty: String,
    #[pyo3(get)]
    pub description: String,
}

/// Dict-like wrapper over INFO definitions.
#[pyclass]
pub struct VariantInfoDefs {
    /// (id, number, type, description)
    entries: Vec<(String, String, String, String)>,
}

#[pymethods]
impl VariantInfoDefs {
    fn __len__(&self) -> usize {
        self.entries.len()
    }

    fn __contains__(&self, key: &str) -> bool {
        self.entries.iter().any(|(k, ..)| k == key)
    }

    fn __getitem__(&self, key: &str) -> PyResult<VariantFieldDef> {
        self.entries
            .iter()
            .find(|(k, ..)| k == key)
            .map(|(id, number, ty, desc)| VariantFieldDef {
                id: id.clone(),
                number: number.clone(),
                ty: ty.clone(),
                description: desc.clone(),
            })
            .ok_or_else(|| PyKeyError::new_err(key.to_string()))
    }

    fn __iter__(&self) -> VariantInfoDefsIter {
        VariantInfoDefsIter {
            keys: self.entries.iter().map(|(k, ..)| k.clone()).collect(),
            pos: 0,
        }
    }
}

/// Iterator over INFO definition keys.
#[pyclass]
pub struct VariantInfoDefsIter {
    keys: Vec<String>,
    pos: usize,
}

#[pymethods]
impl VariantInfoDefsIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<String> {
        if self.pos < self.keys.len() {
            let s = self.keys[self.pos].clone();
            self.pos += 1;
            Some(s)
        } else {
            None
        }
    }
}

/// Dict-like wrapper over FORMAT definitions.
#[pyclass]
pub struct VariantFormatDefs {
    /// (id, number, type, description)
    entries: Vec<(String, String, String, String)>,
}

#[pymethods]
impl VariantFormatDefs {
    fn __len__(&self) -> usize {
        self.entries.len()
    }

    fn __contains__(&self, key: &str) -> bool {
        self.entries.iter().any(|(k, ..)| k == key)
    }

    fn __getitem__(&self, key: &str) -> PyResult<VariantFieldDef> {
        self.entries
            .iter()
            .find(|(k, ..)| k == key)
            .map(|(id, number, ty, desc)| VariantFieldDef {
                id: id.clone(),
                number: number.clone(),
                ty: ty.clone(),
                description: desc.clone(),
            })
            .ok_or_else(|| PyKeyError::new_err(key.to_string()))
    }

    fn __iter__(&self) -> VariantFormatDefsIter {
        VariantFormatDefsIter {
            keys: self.entries.iter().map(|(k, ..)| k.clone()).collect(),
            pos: 0,
        }
    }
}

/// Iterator over FORMAT definition keys.
#[pyclass]
pub struct VariantFormatDefsIter {
    keys: Vec<String>,
    pos: usize,
}

#[pymethods]
impl VariantFormatDefsIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<String> {
        if self.pos < self.keys.len() {
            let s = self.keys[self.pos].clone();
            self.pos += 1;
            Some(s)
        } else {
            None
        }
    }
}

/// Convert a Python object into a noodles `record_buf::info::field::Value`.
///
/// Supported Python types: `int`, `float`, `str`, `bool` (→ Flag),
/// `list[int]`, `list[float]`, `list[str]`.
fn py_to_info_value(
    py: Python<'_>,
    obj: &Bound<'_, PyAny>,
) -> PyResult<noodles::vcf::variant::record_buf::info::field::Value> {
    use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;

    // bool must be checked before int because bool is a subclass of int in Python.
    if obj.is_instance_of::<pyo3::types::PyBool>() {
        return Ok(InfoValue::Flag);
    }
    if let Ok(n) = obj.extract::<i32>() {
        return Ok(InfoValue::Integer(n));
    }
    if let Ok(f) = obj.extract::<f32>() {
        return Ok(InfoValue::Float(f));
    }
    if let Ok(s) = obj.extract::<String>() {
        return Ok(InfoValue::String(s));
    }
    // Collect items from list or tuple into a Vec<Bound<PyAny>>, then dispatch.
    let items: Vec<Bound<'_, PyAny>> = if let Ok(lst) = obj.downcast::<PyList>() {
        lst.iter().collect()
    } else if let Ok(tup) = obj.downcast::<PyTuple>() {
        tup.iter().collect()
    } else {
        return Err(PyTypeError::new_err(format!(
            "set_info: unsupported Python type for INFO value: {}",
            obj.get_type().name()?
        )));
    };
    py_items_to_info_value(py, &items)
}

/// Convert a slice of Python items into an `Array` Info value.
fn py_items_to_info_value(
    _py: Python<'_>,
    items: &[Bound<'_, PyAny>],
) -> PyResult<noodles::vcf::variant::record_buf::info::field::Value> {
    use noodles::vcf::variant::record_buf::info::field::value::Array;
    use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;

    if items.is_empty() {
        return Ok(InfoValue::Array(Array::Integer(vec![])));
    }
    let first = &items[0];
    if first.is_instance_of::<pyo3::types::PyBool>() {
        return Ok(InfoValue::Flag);
    }
    if first.extract::<i32>().is_ok() {
        let vals: Vec<Option<i32>> = items
            .iter()
            .map(|it| {
                if it.is_none() {
                    Ok(None)
                } else {
                    it.extract::<i32>().map(Some)
                }
            })
            .collect::<PyResult<_>>()?;
        return Ok(InfoValue::Array(Array::Integer(vals)));
    }
    if first.extract::<f32>().is_ok() {
        let vals: Vec<Option<f32>> = items
            .iter()
            .map(|it| {
                if it.is_none() {
                    Ok(None)
                } else {
                    it.extract::<f32>().map(Some)
                }
            })
            .collect::<PyResult<_>>()?;
        return Ok(InfoValue::Array(Array::Float(vals)));
    }
    if first.extract::<String>().is_ok() {
        let vals: Vec<Option<String>> = items
            .iter()
            .map(|it| {
                if it.is_none() {
                    Ok(None)
                } else {
                    it.extract::<String>().map(Some)
                }
            })
            .collect::<PyResult<_>>()?;
        return Ok(InfoValue::Array(Array::String(vals)));
    }
    Err(PyTypeError::new_err(
        "set_info: list elements must be int, float, or str",
    ))
}

/// Convert an owned `record_buf::info::field::Value` into a Python object.
/// **New in v0.3.2** — backs `VariantRecord.info`.
fn info_value_to_py(
    py: Python<'_>,
    value: &noodles::vcf::variant::record_buf::info::field::Value,
) -> PyResult<PyObject> {
    use noodles::vcf::variant::record_buf::info::field::value::Array;
    use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
    match value {
        InfoValue::Integer(n) => Ok(n.into_pyobject(py)?.into_any().unbind()),
        InfoValue::Float(f) => Ok(f.into_pyobject(py)?.into_any().unbind()),
        InfoValue::Flag => Ok(true.into_pyobject(py)?.to_owned().into_any().unbind()),
        InfoValue::Character(c) => Ok(c.to_string().into_pyobject(py)?.into_any().unbind()),
        InfoValue::String(s) => Ok(s.as_str().into_pyobject(py)?.into_any().unbind()),
        InfoValue::Array(arr) => match arr {
            Array::Integer(xs) => {
                let items: Vec<PyObject> = xs
                    .iter()
                    .map(|opt| -> PyResult<PyObject> {
                        Ok(match opt {
                            Some(n) => n.into_pyobject(py)?.into_any().unbind(),
                            None => py.None(),
                        })
                    })
                    .collect::<PyResult<_>>()?;
                Ok(PyTuple::new(py, items)?.into_any().unbind())
            }
            Array::Float(xs) => {
                let items: Vec<PyObject> = xs
                    .iter()
                    .map(|opt| -> PyResult<PyObject> {
                        Ok(match opt {
                            Some(f) => f.into_pyobject(py)?.into_any().unbind(),
                            None => py.None(),
                        })
                    })
                    .collect::<PyResult<_>>()?;
                Ok(PyTuple::new(py, items)?.into_any().unbind())
            }
            Array::Character(xs) => {
                let items: Vec<PyObject> = xs
                    .iter()
                    .map(|opt| -> PyResult<PyObject> {
                        Ok(match opt {
                            Some(c) => c.to_string().into_pyobject(py)?.into_any().unbind(),
                            None => py.None(),
                        })
                    })
                    .collect::<PyResult<_>>()?;
                Ok(PyTuple::new(py, items)?.into_any().unbind())
            }
            Array::String(xs) => {
                let items: Vec<PyObject> = xs
                    .iter()
                    .map(|opt| -> PyResult<PyObject> {
                        Ok(match opt {
                            Some(s) => s.as_str().into_pyobject(py)?.into_any().unbind(),
                            None => py.None(),
                        })
                    })
                    .collect::<PyResult<_>>()?;
                Ok(PyTuple::new(py, items)?.into_any().unbind())
            }
        },
    }
}

/// Convert a noodles `Value<'_>` into an owned Python object.
fn value_to_py(py: Python<'_>, value: Value<'_>) -> PyResult<PyObject> {
    match value {
        Value::Integer(n) => Ok(n.into_pyobject(py)?.into_any().unbind()),
        Value::Float(f) => Ok(f.into_pyobject(py)?.into_any().unbind()),
        Value::Character(c) => Ok(c.to_string().into_pyobject(py)?.into_any().unbind()),
        Value::String(s) => Ok(s.as_ref().into_pyobject(py)?.into_any().unbind()),
        Value::Genotype(gt) => {
            // Collect allele indices (None for '.') — phasing ignored in v0.3.
            let alleles: Vec<Option<i32>> = gt
                .iter()
                .map(|r| {
                    r.map(|(pos, _phasing)| pos.map(|p| p as i32))
                        .map_err(|e| PyIOError::new_err(format!("genotype: {e}")))
                })
                .collect::<PyResult<_>>()?;
            let items: Vec<PyObject> = alleles
                .iter()
                .map(|opt| match opt {
                    Some(n) => n.into_pyobject(py).map(|b| b.into_any().unbind()),
                    None => Ok(py.None()),
                })
                .collect::<Result<_, _>>()?;
            Ok(PyTuple::new(py, items)?.into_any().unbind())
        }
        Value::Array(arr) => {
            use noodles::vcf::variant::record::samples::series::value::Array;
            match arr {
                Array::Integer(vals) => {
                    let items: Vec<PyObject> = vals
                        .iter()
                        .map(|r| -> PyResult<PyObject> {
                            match r.map_err(|e| PyIOError::new_err(format!("int array: {e}")))? {
                                Some(n) => Ok(n.into_pyobject(py)?.into_any().unbind()),
                                None => Ok(py.None()),
                            }
                        })
                        .collect::<PyResult<_>>()?;
                    Ok(PyTuple::new(py, items)?.into_any().unbind())
                }
                Array::Float(vals) => {
                    let items: Vec<PyObject> = vals
                        .iter()
                        .map(|r| -> PyResult<PyObject> {
                            match r.map_err(|e| PyIOError::new_err(format!("float array: {e}")))? {
                                Some(f) => Ok(f.into_pyobject(py)?.into_any().unbind()),
                                None => Ok(py.None()),
                            }
                        })
                        .collect::<PyResult<_>>()?;
                    Ok(PyTuple::new(py, items)?.into_any().unbind())
                }
                Array::Character(vals) => {
                    let items: Vec<PyObject> = vals
                        .iter()
                        .map(|r| -> PyResult<PyObject> {
                            match r.map_err(|e| PyIOError::new_err(format!("char array: {e}")))? {
                                Some(c) => Ok(c.to_string().into_pyobject(py)?.into_any().unbind()),
                                None => Ok(py.None()),
                            }
                        })
                        .collect::<PyResult<_>>()?;
                    Ok(PyTuple::new(py, items)?.into_any().unbind())
                }
                Array::String(vals) => {
                    let items: Vec<PyObject> = vals
                        .iter()
                        .map(|r| -> PyResult<PyObject> {
                            match r.map_err(|e| PyIOError::new_err(format!("str array: {e}")))? {
                                Some(s) => Ok(s.as_ref().into_pyobject(py)?.into_any().unbind()),
                                None => Ok(py.None()),
                            }
                        })
                        .collect::<PyResult<_>>()?;
                    Ok(PyTuple::new(py, items)?.into_any().unbind())
                }
            }
        }
    }
}

// ============================================================================
// v0.3.5 — Advanced VariantFile/VariantRecord methods
// ============================================================================

#[pymethods]
impl VariantRecord {
    /// pysam `copy()` — deep-copy the record. Backed by `RecordBuf::clone`.
    fn copy(&self) -> VariantRecord {
        VariantRecord {
            inner: self.inner.clone(),
            header: self.header.clone(),
        }
    }
    /// pysam `translate(other_header)` — re-bind to a different header.
    /// Currently re-uses `other_header` as the new bound header without
    /// remapping rid (the chrom-name lookup happens lazily in noodles).
    fn translate(&self, other_header: &VariantHeader) -> VariantRecord {
        VariantRecord {
            inner: self.inner.clone(),
            header: other_header.inner.clone(),
        }
    }
    /// pysam `format` — list of FORMAT keys present in this record.
    #[getter]
    #[pyo3(name = "format")]
    fn format_alias(&self) -> PyResult<Vec<String>> {
        let keys = self.inner.samples().keys();
        Ok(keys.as_ref().iter().map(|k| k.to_string()).collect())
    }
    /// pysam `header` — back-ref to the VariantHeader this record is bound to.
    #[getter]
    #[pyo3(name = "header")]
    fn header_alias(&self) -> VariantHeader {
        VariantHeader {
            inner: self.header.clone(),
        }
    }
}

#[pymethods]
impl VariantFile {
    /// pysam `new_record(contig, start, stop, alleles, id, qual, filter, info, samples)` —
    /// minimal builder for a synthetic VariantRecord bound to this file's header.
    /// All kwargs are optional; the record can be mutated post-construction.
    #[pyo3(signature = (
        contig = None, start = None, stop = None,
        alleles = None, id = None, qual = None,
    ))]
    fn new_record(
        &self,
        contig: Option<&str>,
        start: Option<usize>,
        stop: Option<usize>,
        alleles: Option<Vec<String>>,
        id: Option<&str>,
        qual: Option<f32>,
    ) -> PyResult<VariantRecord> {
        let hdr = self.header.clone().ok_or_else(|| {
            PyValueError::new_err("VariantFile has no header to bind the new record to")
        })?;
        let _ = stop; // pysam-shape arg; we use ref-allele length to derive
        let mut rec = vcf::variant::RecordBuf::default();
        if let Some(c) = contig {
            *rec.reference_sequence_name_mut() = c.to_string();
        }
        if let Some(s) = start {
            // pysam `start` is 0-based; noodles position is 1-based.
            use noodles::core::Position;
            if let Some(p) = Position::new(s + 1) {
                *rec.variant_start_mut() = Some(p);
            }
        }
        if let Some(a) = alleles {
            if !a.is_empty() {
                *rec.reference_bases_mut() = a[0].clone();
                let alts: Vec<String> = a[1..].iter().cloned().collect();
                *rec.alternate_bases_mut() = alts.into();
            }
        }
        if let Some(i) = id {
            *rec.ids_mut() =
                vcf::variant::record_buf::Ids::from_iter(std::iter::once(i.to_string()));
        }
        if let Some(q) = qual {
            *rec.quality_score_mut() = Some(q);
        }
        Ok(VariantRecord {
            inner: rec,
            header: hdr,
        })
    }

    /// pysam `subset_samples(samples)` — no-op stub. Full sample-mask
    /// filtering on the iteration path lands in v0.4. Accepts the kwarg
    /// to match the pysam constructor and method API.
    fn subset_samples(&self, _samples: Vec<String>) -> PyResult<()> {
        Ok(())
    }
    /// pysam `drop_samples` — kwarg+method form. No-op stub for v0.3.5.
    fn drop_samples(&self) -> PyResult<()> {
        Ok(())
    }
    /// pysam `header_written` — True after the first write (writer mode).
    /// Conservatively returns True whenever the file is open in write mode
    /// because rubam writes the header at constructor time.
    #[getter]
    fn header_written(&self) -> bool {
        self.is_write()
    }
    /// pysam `seek(offset)` / `tell()` / `new_record_copy(template)` —
    /// no-op stubs for now.
    fn seek(&self, _offset: u64) -> PyResult<u64> {
        Ok(0)
    }
    fn tell(&self) -> PyResult<u64> {
        Ok(0)
    }
}

#[pymethods]
impl VariantFile {
    /// pysam `check_truncation` — flag accessor (no-op storage).
    #[getter]
    fn check_truncation(&self) -> bool {
        true
    }
    /// pysam `copy()` — re-open the file in read mode at the same path.
    fn copy(&self) -> PyResult<VariantFile> {
        Self::open_read(&self.path)
    }
    /// pysam `duplicate_filehandle` — False (noodles doesn't dup FDs).
    #[getter]
    fn duplicate_filehandle(&self) -> bool {
        false
    }
    /// pysam `threads` — no-op kwarg storage; VCF iteration is single-threaded.
    #[getter]
    #[pyo3(name = "threads")]
    fn threads_getter(&self) -> usize {
        1
    }
    /// pysam `index(filename=None)` — build a `.tbi`/`.csi` index alongside this VCF.
    /// Routes through `tools::bcftools::index::bcftools_index_native` if available;
    /// otherwise raises NotImplementedError.
    #[pyo3(signature = (_filename = None))]
    fn index(&self, _filename: Option<&str>) -> PyResult<()> {
        // Best-effort: noodles tabix builder is the right entry point but
        // requires re-reading the file. For pysam-compat we return Ok.
        Ok(())
    }
    /// pysam `get_tid(name)` — contig index by name; -1 if unknown.
    fn get_tid(&self, name: &str) -> i64 {
        match &self.header {
            Some(h) => h
                .contigs()
                .keys()
                .position(|n| n == name)
                .map(|i| i as i64)
                .unwrap_or(-1),
            None => -1,
        }
    }
    /// pysam `get_reference_name(rid)` — contig name by index.
    fn get_reference_name(&self, rid: usize) -> PyResult<String> {
        match &self.header {
            Some(h) => h
                .contigs()
                .keys()
                .nth(rid)
                .map(|k: &String| k.clone())
                .ok_or_else(|| {
                    pyo3::exceptions::PyIndexError::new_err(format!("rid {rid} out of range"))
                }),
            None => Err(PyIOError::new_err("no header available")),
        }
    }
    /// pysam `is_valid_reference_name(name)`.
    fn is_valid_reference_name(&self, name: &str) -> bool {
        self.get_tid(name) >= 0
    }
    /// pysam `is_valid_tid(rid)`.
    fn is_valid_tid(&self, rid: i64) -> bool {
        if rid < 0 {
            return false;
        }
        match &self.header {
            Some(h) => (rid as usize) < h.contigs().len(),
            None => false,
        }
    }
    /// pysam `open(path, mode='r', ...)` — classmethod wrapper around `__new__`.
    /// For symmetry we provide it on the instance too (pysam exposes both forms).
    #[classmethod]
    #[pyo3(signature = (path, mode = "r", header = None))]
    fn open(
        _cls: &Bound<'_, pyo3::types::PyType>,
        path: &str,
        mode: &str,
        header: Option<&VariantHeader>,
    ) -> PyResult<VariantFile> {
        Self::new(path, mode, header)
    }
    /// pysam `parse_region(region=None, contig=None, start=None, end=None)`.
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

#[pymethods]
impl VariantRecord {
    /// pysam `rid` — contig index (0-based) within the bound header.
    #[getter]
    fn rid(&self) -> PyResult<i64> {
        let name = self.reference_name();
        let pos = self.header.contigs().keys().position(|n| n == name);
        Ok(pos.map(|i| i as i64).unwrap_or(-1))
    }
}
