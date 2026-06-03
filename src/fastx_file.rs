//! `rubam.FastxFile` — pysam-compatible streaming FASTA/FASTQ reader.
//!
//! Auto-detects:
//!   - `.fa` / `.fasta` / `.fna`  → pure FASTA (no quality)
//!   - `.fq` / `.fastq`           → FASTQ
//!   - `.fa.gz` / `.fastq.gz`     → BGZF or gzip wrapped
//!
//! ```python
//! import rubam
//! with rubam.FastxFile("reads.fastq.gz") as fx:
//!     for entry in fx:
//!         print(entry.name, entry.sequence, entry.quality, entry.comment)
//! ```

use std::cell::RefCell;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::PathBuf;

use noodles::fastq;

use pyo3::exceptions::{PyIOError, PyStopIteration};
use pyo3::prelude::*;

#[pyclass(unsendable)]
pub struct FastxRecord {
    name: String,
    sequence: String,
    quality: Option<String>,
    comment: Option<String>,
}

#[pymethods]
impl FastxRecord {
    #[getter]
    fn name(&self) -> String {
        self.name.clone()
    }
    #[getter]
    fn sequence(&self) -> String {
        self.sequence.clone()
    }
    /// pysam `quality` — phred+33-encoded string; None for FASTA.
    #[getter]
    fn quality(&self) -> Option<String> {
        self.quality.clone()
    }
    #[getter]
    fn comment(&self) -> Option<String> {
        self.comment.clone()
    }
    fn __repr__(&self) -> String {
        format!(
            "FastxRecord(name={:?}, len={}, has_quality={})",
            self.name,
            self.sequence.len(),
            self.quality.is_some()
        )
    }
}

enum FastxReader {
    /// FASTA: we drive a plain BufReader and stash the next header line
    /// so that records terminate cleanly on the next `>` boundary.
    Fasta {
        reader: BufReader<Box<dyn Read + Send>>,
        next_header: Option<String>,
    },
    Fastq(fastq::io::Reader<BufReader<Box<dyn Read + Send>>>),
}

#[pyclass(unsendable)]
pub struct FastxFile {
    path: PathBuf,
    inner: RefCell<Option<FastxReader>>,
}

fn open_reader(path: &PathBuf) -> PyResult<FastxReader> {
    let lc = path.to_string_lossy().to_ascii_lowercase();
    let f = File::open(path)
        .map_err(|e| PyIOError::new_err(format!("failed to open {}: {e}", path.display())))?;
    let inner: Box<dyn Read + Send> = if lc.ends_with(".gz") {
        Box::new(flate2::read::MultiGzDecoder::new(f))
    } else {
        Box::new(f)
    };
    let buf = BufReader::new(inner);
    if lc.contains(".fastq") || lc.contains(".fq") {
        Ok(FastxReader::Fastq(fastq::io::Reader::new(buf)))
    } else {
        Ok(FastxReader::Fasta {
            reader: buf,
            next_header: None,
        })
    }
}

#[pymethods]
impl FastxFile {
    #[new]
    fn new(path: PathBuf) -> PyResult<Self> {
        let r = open_reader(&path)?;
        Ok(FastxFile {
            path,
            inner: RefCell::new(Some(r)),
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&self) -> PyResult<FastxRecord> {
        let mut inner_opt = self.inner.borrow_mut();
        let r = inner_opt
            .as_mut()
            .ok_or_else(|| PyIOError::new_err("FastxFile is closed"))?;
        match r {
            FastxReader::Fasta {
                reader,
                next_header,
            } => {
                // Start with either a stashed header or the next `>` line.
                let header_line = if let Some(h) = next_header.take() {
                    h
                } else {
                    let mut line = String::new();
                    loop {
                        line.clear();
                        let n = reader
                            .read_line(&mut line)
                            .map_err(|e| PyIOError::new_err(format!("fasta read: {e}")))?;
                        if n == 0 {
                            return Err(PyStopIteration::new_err(""));
                        }
                        let t = line.trim_end();
                        if t.is_empty() {
                            continue;
                        }
                        if t.starts_with('>') {
                            break;
                        }
                        // We are between records — silently skip non-blank
                        // non-header lines until the next `>`.
                    }
                    line.clone()
                };
                let mut seq = String::new();
                loop {
                    let mut line = String::new();
                    let n = reader
                        .read_line(&mut line)
                        .map_err(|e| PyIOError::new_err(format!("fasta read: {e}")))?;
                    if n == 0 {
                        break;
                    }
                    let t = line.trim_end();
                    if t.starts_with('>') {
                        // Stash for the next call.
                        *next_header = Some(line);
                        break;
                    }
                    seq.push_str(t);
                }
                let header = header_line.trim_start_matches('>').trim_end().to_string();
                let (name, comment) = match header.split_once(char::is_whitespace) {
                    Some((n, c)) => (n.to_string(), Some(c.to_string())),
                    None => (header, None),
                };
                Ok(FastxRecord {
                    name,
                    sequence: seq,
                    quality: None,
                    comment,
                })
            }
            FastxReader::Fastq(rdr) => {
                let mut rec = fastq::Record::default();
                let n = rdr
                    .read_record(&mut rec)
                    .map_err(|e| PyIOError::new_err(format!("fastq read: {e}")))?;
                if n == 0 {
                    return Err(PyStopIteration::new_err(""));
                }
                let name = String::from_utf8_lossy(rec.name()).into_owned();
                let desc = rec.description();
                let comment = if desc.is_empty() {
                    None
                } else {
                    Some(String::from_utf8_lossy(desc).into_owned())
                };
                let sequence = String::from_utf8_lossy(rec.sequence()).into_owned();
                let quality = Some(String::from_utf8_lossy(rec.quality_scores()).into_owned());
                Ok(FastxRecord {
                    name,
                    sequence,
                    quality,
                    comment,
                })
            }
        }
    }

    fn close(&self) -> PyResult<()> {
        *self.inner.borrow_mut() = None;
        Ok(())
    }
    #[getter]
    fn is_open(&self) -> bool {
        self.inner.borrow().is_some()
    }
    #[getter]
    fn filename(&self) -> String {
        self.path.to_string_lossy().into_owned()
    }

    fn __enter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __exit__(
        &self,
        _exc_type: Option<PyObject>,
        _exc_value: Option<PyObject>,
        _exc_tb: Option<PyObject>,
    ) -> PyResult<bool> {
        self.close()?;
        Ok(false)
    }
}
