//! `TabixFile` — pysam-compatible tabix-indexed file reader (v0.3.3 stub).
//!
//! v0.3.3 ships the *class shape* so pysam-porting code that does
//! `isinstance(x, rubam.TabixFile)` doesn't AttributeError, but the
//! constructor raises `NotImplementedError`. A full implementation is
//! tracked for v0.4 once the noodles 0.107 tabix reader API stabilises
//! (the existing prototype hit several private-trait import paths and
//! BGZF virtual-position helpers that the 0.107 release does not
//! expose publicly).
//!
//! Until then, the supported alternatives are:
//! - `rubam.VariantFile(path)` for `.vcf.gz` / `.bcf` random access (uses tbi/csi internally),
//! - subprocess to system `tabix` for arbitrary tabix-indexed files.

#![cfg(feature = "python")]

use pyo3::exceptions::PyNotImplementedError;
use pyo3::prelude::*;

#[pyclass(module = "rubam", unsendable)]
pub struct TabixFile;

#[pymethods]
impl TabixFile {
    #[new]
    #[pyo3(signature = (_filename, *, _index = None))]
    fn new(_filename: &str, _index: Option<&str>) -> PyResult<Self> {
        Err(PyNotImplementedError::new_err(
            "rubam.TabixFile is reserved for a v0.4 implementation. \
             Until then, use rubam.VariantFile for .vcf.gz / .bcf \
             random access, or shell out to `tabix` via subprocess. \
             See docs/pysam_compatibility_matrix.md.",
        ))
    }
}

#[pyclass(module = "rubam", unsendable)]
pub struct TabixFileIter;

#[pymethods]
impl TabixFileIter {
    fn __iter__(slf: PyRef<Self>) -> PyRef<Self> {
        slf
    }

    fn __next__(&self) -> PyResult<Option<String>> {
        Err(PyNotImplementedError::new_err(
            "rubam.TabixFile is reserved for v0.4",
        ))
    }
}
