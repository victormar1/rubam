//! samtools faidx — build the .fai index and optionally extract a subsequence.

use std::fs::File;
use std::io::BufWriter;

use noodles::core::Region;
use noodles::fasta;
#[cfg(feature = "python")]
use pyo3::exceptions::{PyIOError, PyValueError};
#[cfg(feature = "python")]
use pyo3::prelude::*;

/// `samtools faidx` equivalent. Builds `<input>.fai` if missing.
///
/// If `region` is given (`chr:start-end`, 1-based inclusive), returns the
/// subsequence as an upper-case string. Otherwise returns None and only
/// writes the index file.
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (input, *, region = None))]
pub fn faidx(input: &str, region: Option<&str>) -> PyResult<Option<String>> {
    faidx_index_only(input).map_err(|e| PyIOError::new_err(format!("faidx index: {e}")))?;
    let Some(r) = region else {
        return Ok(None);
    };
    let (header, seq) = faidx_subseq(input, r).map_err(|e| match e.kind() {
        std::io::ErrorKind::InvalidInput => PyValueError::new_err(e.to_string()),
        _ => PyIOError::new_err(format!("faidx: {e}")),
    })?;
    let _ = header;
    Ok(Some(seq))
}

/// Build the .fai if it doesn't exist; no-op otherwise.
pub fn faidx_index_only(input: &str) -> std::io::Result<()> {
    let fai_path = format!("{}.fai", input);
    if std::path::Path::new(&fai_path).exists() {
        return Ok(());
    }
    let index = fasta::fs::index(input)?;
    let f = File::create(&fai_path)?;
    let mut writer = fasta::fai::io::Writer::new(BufWriter::new(f));
    writer.write_index(&index)?;
    Ok(())
}

/// Extract a subsequence given a `chr:start-end` 1-based-inclusive region.
/// Returns (display_header, sequence). `display_header` matches `samtools faidx`
/// output style: `chr:start-end`.
pub fn faidx_subseq(input: &str, region_str: &str) -> std::io::Result<(String, String)> {
    // Validate the region string before building the reader
    let region: Region = region_str.parse().map_err(|_| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!("invalid region {region_str:?}: must be 'chrom:start-end'"),
        )
    })?;
    faidx_index_only(input)?;
    let mut reader = fasta::io::indexed_reader::Builder::default().build_from_path(input)?;
    let record = reader.query(&region)?;
    let bytes = record.sequence().as_ref();
    let seq = String::from_utf8_lossy(bytes).into_owned();
    Ok((region_str.to_string(), seq))
}
