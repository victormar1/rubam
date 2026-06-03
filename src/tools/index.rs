//! samtools index — write a .bai (BAI) for a coordinate-sorted BAM.

use std::fs::File;

use noodles::bam;
#[cfg(feature = "python")]
use pyo3::exceptions::PyIOError;
#[cfg(feature = "python")]
use pyo3::prelude::*;

/// Build the BAI index for a coordinate-sorted BAM. Equivalent to
/// `samtools index INPUT`, writing alongside as `INPUT.bai`.
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (input, *, csi = false))]
pub fn index(input: &str, csi: bool) -> PyResult<()> {
    if csi {
        return Err(PyIOError::new_err("CSI indexing lands in v0.2.x"));
    }
    index_native(input).map_err(|e| PyIOError::new_err(format!("index: {e}")))
}

/// Pure-Rust entry point used by the rubam-samtools shadow CLI in Phase C.
pub fn index_native(input: &str) -> std::io::Result<()> {
    let index = bam::fs::index(input)?;
    let bai_path = format!("{}.bai", input);
    let bai_file = File::create(&bai_path)?;
    let mut writer = bam::bai::io::Writer::new(bai_file);
    writer.write_index(&index)?;
    Ok(())
}
