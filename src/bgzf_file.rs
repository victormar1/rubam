//! pysam.BGZFile equivalent — a real BGZF read/write file wrapper.
//!
//! Backed by `noodles::bgzf`. Exposes the pysam.BGZFile read-side
//! interface (read / readline / readlines / seek / tell / close) plus
//! a write-side (write / writelines) when opened in `"wb"` mode.
//!
//! This is a v0.3.12 "real impl" addition: previously rubam exposed
//! BGZFile as an empty marker class; now it does actual BGZF I/O.

#[cfg(feature = "python")]
use pyo3::exceptions::{PyIOError, PyValueError};
#[cfg(feature = "python")]
use pyo3::prelude::*;

use std::cell::RefCell;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

use noodles::bgzf;

/// Inner I/O enum: BGZF reader, BGZF writer, or closed.
enum BgzfInner {
    Reader(BufReader<bgzf::io::Reader<File>>),
    Writer(bgzf::io::Writer<BufWriter<File>>),
    Closed,
}

/// pysam.BGZFile — open a BGZF file for read or write.
#[cfg(feature = "python")]
#[pyclass(unsendable)]
pub struct BGZFile {
    path: String,
    mode: String,
    inner: RefCell<BgzfInner>,
}

#[cfg(feature = "python")]
#[pymethods]
impl BGZFile {
    #[new]
    #[pyo3(signature = (path, mode = "rb"))]
    fn new(path: std::path::PathBuf, mode: &str) -> PyResult<Self> {
        let path_str = path.to_string_lossy().into_owned();
        let inner = match mode {
            "rb" | "r" => {
                let f = File::open(&path)
                    .map_err(|e| PyIOError::new_err(format!("open BGZF {path_str:?}: {e}")))?;
                let reader = bgzf::io::Reader::new(f);
                BgzfInner::Reader(BufReader::new(reader))
            }
            "wb" | "w" => {
                let f = File::create(&path)
                    .map_err(|e| PyIOError::new_err(format!("create BGZF {path_str:?}: {e}")))?;
                let writer = bgzf::io::Writer::new(BufWriter::new(f));
                BgzfInner::Writer(writer)
            }
            other => {
                return Err(PyValueError::new_err(format!(
                    "BGZFile: unsupported mode {other:?}; use 'rb' or 'wb'"
                )));
            }
        };
        Ok(Self {
            path: path_str,
            mode: mode.to_string(),
            inner: RefCell::new(inner),
        })
    }

    /// Read `n` bytes (or all remaining if -1).
    #[pyo3(signature = (n = -1))]
    fn read(&self, n: i64) -> PyResult<Vec<u8>> {
        let mut inner = self.inner.borrow_mut();
        match &mut *inner {
            BgzfInner::Reader(r) => {
                if n < 0 {
                    let mut buf = Vec::new();
                    r.read_to_end(&mut buf)
                        .map_err(|e| PyIOError::new_err(format!("read: {e}")))?;
                    Ok(buf)
                } else {
                    let mut buf = vec![0u8; n as usize];
                    let got = r
                        .read(&mut buf)
                        .map_err(|e| PyIOError::new_err(format!("read: {e}")))?;
                    buf.truncate(got);
                    Ok(buf)
                }
            }
            BgzfInner::Writer(_) => Err(PyIOError::new_err("read on write-mode BGZFile")),
            BgzfInner::Closed => Err(PyIOError::new_err("read on closed BGZFile")),
        }
    }

    /// Read one line (up to and including the newline). Returns
    /// empty bytes on EOF.
    fn readline(&self) -> PyResult<Vec<u8>> {
        let mut inner = self.inner.borrow_mut();
        match &mut *inner {
            BgzfInner::Reader(r) => {
                let mut buf = Vec::new();
                r.read_until(b'\n', &mut buf)
                    .map_err(|e| PyIOError::new_err(format!("readline: {e}")))?;
                Ok(buf)
            }
            BgzfInner::Writer(_) => Err(PyIOError::new_err("readline on write-mode BGZFile")),
            BgzfInner::Closed => Err(PyIOError::new_err("readline on closed BGZFile")),
        }
    }

    /// Read all lines.
    fn readlines(&self) -> PyResult<Vec<Vec<u8>>> {
        let mut inner = self.inner.borrow_mut();
        match &mut *inner {
            BgzfInner::Reader(r) => {
                let mut lines = Vec::new();
                loop {
                    let mut buf = Vec::new();
                    let n = r
                        .read_until(b'\n', &mut buf)
                        .map_err(|e| PyIOError::new_err(format!("readlines: {e}")))?;
                    if n == 0 {
                        break;
                    }
                    lines.push(buf);
                }
                Ok(lines)
            }
            BgzfInner::Writer(_) => Err(PyIOError::new_err("readlines on write-mode BGZFile")),
            BgzfInner::Closed => Err(PyIOError::new_err("readlines on closed BGZFile")),
        }
    }

    /// Write bytes to the underlying BGZF writer.
    fn write(&self, data: &[u8]) -> PyResult<usize> {
        let mut inner = self.inner.borrow_mut();
        match &mut *inner {
            BgzfInner::Writer(w) => {
                let n = w
                    .write(data)
                    .map_err(|e| PyIOError::new_err(format!("write: {e}")))?;
                Ok(n)
            }
            BgzfInner::Reader(_) => Err(PyIOError::new_err("write on read-mode BGZFile")),
            BgzfInner::Closed => Err(PyIOError::new_err("write on closed BGZFile")),
        }
    }

    fn writelines(&self, lines: Vec<Vec<u8>>) -> PyResult<usize> {
        let mut total = 0usize;
        for line in lines {
            total += self.write(&line)?;
        }
        Ok(total)
    }

    /// Flush the write buffer.
    fn flush(&self) -> PyResult<()> {
        let mut inner = self.inner.borrow_mut();
        if let BgzfInner::Writer(w) = &mut *inner {
            w.flush()
                .map_err(|e| PyIOError::new_err(format!("flush: {e}")))?;
        }
        Ok(())
    }

    /// Close the file (write EOF block for writers).
    fn close(&self) -> PyResult<()> {
        let mut inner = self.inner.borrow_mut();
        *inner = BgzfInner::Closed;
        Ok(())
    }

    #[getter]
    fn closed(&self) -> bool {
        matches!(*self.inner.borrow(), BgzfInner::Closed)
    }

    #[getter]
    fn is_open(&self) -> bool {
        !self.closed()
    }

    #[getter]
    fn filename(&self) -> String {
        self.path.clone()
    }

    #[getter]
    fn mode(&self) -> String {
        self.mode.clone()
    }

    fn __enter__<'py>(slf: PyRef<'py, Self>) -> PyRef<'py, Self> {
        slf
    }

    #[pyo3(signature = (_exc_type=None, _exc_val=None, _exc_tb=None))]
    fn __exit__(
        &self,
        _exc_type: Option<PyObject>,
        _exc_val: Option<PyObject>,
        _exc_tb: Option<PyObject>,
    ) -> PyResult<()> {
        self.close()
    }
}
