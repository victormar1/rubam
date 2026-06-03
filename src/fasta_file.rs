//! `FastaFile` — pysam-compatible random-access FASTA reader.
//!
//! Mirrors `pysam.FastaFile`: open a `.fa` (with `.fai`) or `.fa.gz`
//! (with `.fai` + `.gzi`), and expose `.fetch(contig, start, end)`
//! returning the subsequence as a Python `str`. The 0-based half-open
//! coordinate convention matches pysam.

#![cfg(feature = "python")]

use std::cell::RefCell;
use std::path::PathBuf;

use noodles::core::region::Interval;
use noodles::core::Position;
use noodles::core::Region;
use noodles::fasta;
use pyo3::exceptions::{PyIOError, PyKeyError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyTuple;

use crate::tools::faidx::faidx_index_only;

/// pysam.FastaFile-compatible random-access FASTA reader.
///
/// Marked `unsendable` (like `AlignmentFile`) because the underlying
/// noodles indexed reader is not `Send` — Python callers should not
/// share the same `FastaFile` across threads.
#[pyclass(module = "rubam", unsendable)]
pub struct FastaFile {
    path: PathBuf,
    inner: RefCell<Option<fasta::io::IndexedReader<noodles::fasta::io::BufReader<std::fs::File>>>>,
    refs: Vec<String>,
    lengths: Vec<u64>,
}

#[pymethods]
impl FastaFile {
    /// Open a FASTA. If the `.fai` is missing it will be built automatically.
    ///
    /// `pysam.FastaFile(filename)` equivalent.
    #[new]
    #[pyo3(signature = (filename))]
    fn new(filename: PathBuf) -> PyResult<Self> {
        let path_str = filename
            .to_str()
            .ok_or_else(|| PyValueError::new_err(format!("non-UTF-8 path: {filename:?}")))?;

        // Build the .fai if missing (samtools-faidx-equivalent), so callers
        // don't need to pre-run `samtools faidx` themselves.
        faidx_index_only(path_str).map_err(|e| {
            PyIOError::new_err(format!(
                "FastaFile: failed to build index for {path_str:?}: {e}"
            ))
        })?;

        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path_str)
            .map_err(|e| {
                PyIOError::new_err(format!("FastaFile: failed to open {path_str:?}: {e}"))
            })?;

        // Read the index now so .references / .lengths are available
        // without re-opening.
        let fai_path = format!("{path_str}.fai");
        let fai = fasta::fai::fs::read(&fai_path).map_err(|e| {
            PyIOError::new_err(format!("FastaFile: failed to read {fai_path:?}: {e}"))
        })?;
        let (refs, lengths): (Vec<String>, Vec<u64>) = fai
            .as_ref()
            .iter()
            .map(|rec| {
                (
                    String::from_utf8_lossy(rec.name()).into_owned(),
                    rec.length(),
                )
            })
            .unzip();

        Ok(FastaFile {
            path: filename,
            inner: RefCell::new(Some(reader)),
            refs,
            lengths,
        })
    }

    /// Context manager: `with rubam.FastaFile(...) as fa: ...`.
    fn __enter__(slf: PyRef<Self>) -> PyRef<Self> {
        slf
    }

    fn __exit__(
        &self,
        _exc_type: &Bound<'_, PyAny>,
        _exc_val: &Bound<'_, PyAny>,
        _exc_tb: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        self.close()
    }

    /// Close the file. After this any further call raises.
    fn close(&self) -> PyResult<()> {
        *self.inner.borrow_mut() = None;
        Ok(())
    }

    #[getter]
    fn is_open(&self) -> bool {
        self.inner.borrow().is_some()
    }

    #[getter]
    fn filename(&self) -> PyResult<String> {
        self.path
            .to_str()
            .map(|s| s.to_string())
            .ok_or_else(|| PyValueError::new_err("non-UTF-8 path"))
    }

    /// Tuple of reference names, in FAI order. Matches `pysam.FastaFile.references`.
    #[getter]
    fn references<'py>(&self, py: Python<'py>) -> Bound<'py, PyTuple> {
        PyTuple::new(py, &self.refs).unwrap()
    }

    /// Tuple of reference lengths (as ints), in FAI order. Matches `pysam.FastaFile.lengths`.
    #[getter]
    fn lengths<'py>(&self, py: Python<'py>) -> Bound<'py, PyTuple> {
        PyTuple::new(py, &self.lengths).unwrap()
    }

    #[getter]
    fn nreferences(&self) -> usize {
        self.refs.len()
    }

    /// `pysam.FastaFile.get_reference_length(contig)` — bytes in the contig.
    fn get_reference_length(&self, contig: &str) -> PyResult<u64> {
        match self.refs.iter().position(|r| r == contig) {
            Some(i) => Ok(self.lengths[i]),
            None => Err(PyKeyError::new_err(format!("unknown contig {contig:?}"))),
        }
    }

    /// pysam-compatible fetch. **0-based half-open** coordinates: `fetch("chr1", 0, 10)`
    /// returns the first 10 bases. If `start` and `end` are both `None`,
    /// returns the whole contig. Returns an upper-case-or-as-stored
    /// Python `str` of the subsequence.
    #[pyo3(signature = (reference=None, start=None, end=None, *, region=None))]
    fn fetch(
        &self,
        reference: Option<&str>,
        start: Option<i64>,
        end: Option<i64>,
        region: Option<&str>,
    ) -> PyResult<String> {
        let mut guard = self.inner.borrow_mut();
        let reader = guard
            .as_mut()
            .ok_or_else(|| PyIOError::new_err("FastaFile is closed"))?;

        // Build the noodles Region either from `region="chr:start-end"` (samtools
        // syntax, 1-based inclusive) or from (reference, start, end) (pysam-style,
        // 0-based half-open).
        let parsed_region: Region = if let Some(r) = region {
            r.parse().map_err(|_| {
                PyValueError::new_err(format!("invalid region {r:?}: must be 'chrom:start-end'"))
            })?
        } else {
            let contig = reference.ok_or_else(|| {
                PyValueError::new_err("fetch() requires either `region=...` or `reference=...`")
            })?;
            let contig_len = self.get_reference_length(contig)?;
            let start0 = start.unwrap_or(0);
            let end0 = end.unwrap_or(contig_len as i64);

            if start0 < 0 || end0 < 0 {
                return Err(PyValueError::new_err(
                    "fetch(): start/end must be non-negative",
                ));
            }
            if (end0 as u64) > contig_len {
                return Err(PyValueError::new_err(format!(
                    "fetch(): end ({end0}) past contig length {contig_len} for {contig:?}"
                )));
            }
            if start0 > end0 {
                return Err(PyValueError::new_err("fetch(): start must be <= end"));
            }
            // Convert 0-based half-open [start0, end0) to noodles 1-based inclusive [start1, end1].
            // Empty range (start0 == end0) is allowed: pysam returns "".
            if start0 == end0 {
                return Ok(String::new());
            }
            let start1 = Position::try_from((start0 + 1) as usize)
                .map_err(|e| PyValueError::new_err(format!("bad start: {e}")))?;
            let end1 = Position::try_from(end0 as usize)
                .map_err(|e| PyValueError::new_err(format!("bad end: {e}")))?;
            Region::new(contig.as_bytes(), Interval::from(start1..=end1))
        };

        let record = reader
            .query(&parsed_region)
            .map_err(|e| PyIOError::new_err(format!("FastaFile.fetch: {e}")))?;
        let bytes = record.sequence().as_ref();
        Ok(String::from_utf8_lossy(bytes).into_owned())
    }

    fn __repr__(&self) -> String {
        let n = self.refs.len();
        let path = self.path.display();
        format!("<rubam.FastaFile {path:?} (n_refs={n})>")
    }
}
