//! AlignmentFile.pileup() iterator (buffered; per-read PileupRead in v0.2.x).

use pyo3::prelude::*;

#[pyclass]
pub struct PileupColumn {
    #[pyo3(get)]
    pub reference_name: String,
    #[pyo3(get)]
    pub reference_pos: usize, // 0-based
    #[pyo3(get)]
    pub depth: u32,
    #[pyo3(get)]
    pub a: u32,
    #[pyo3(get)]
    pub c: u32,
    #[pyo3(get)]
    pub g: u32,
    #[pyo3(get)]
    pub t: u32,
    #[pyo3(get)]
    pub n: u32,
}

#[pymethods]
impl PileupColumn {
    #[getter]
    fn nsegments(&self) -> u32 {
        self.depth
    }
}

#[pyclass(unsendable)]
pub struct PileupIter {
    pub(crate) cols: std::cell::RefCell<std::vec::IntoIter<PileupColumn>>,
}

#[pymethods]
impl PileupIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __next__(&self) -> Option<PileupColumn> {
        self.cols.borrow_mut().next()
    }
}
