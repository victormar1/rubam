//! samtools merge — concatenate multiple coordinate-sorted BAMs into one.
//!
//! v0.2: simple sequential read + re-sort all records in memory. K-way merge
//! of pre-sorted streams lands in v0.2.x.

use std::fs::File;
use std::io::BufWriter;

use noodles::bam;
use noodles::bgzf;
use noodles::sam::{alignment::io::Write as _, alignment::RecordBuf};
#[cfg(feature = "python")]
use pyo3::exceptions::{PyIOError, PyValueError};
#[cfg(feature = "python")]
use pyo3::prelude::*;

use crate::common::{open_streaming, read_header_streaming};

type BamWriter = bam::io::Writer<bgzf::io::Writer<BufWriter<File>>>;

/// Concatenate multiple coordinate-sorted BAMs into a single output BAM.
/// Equivalent to `samtools merge -f OUTPUT INPUTS...` (force=true by default).
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (inputs, output, *, force = true))]
pub fn merge(inputs: Vec<String>, output: &str, force: bool) -> PyResult<()> {
    if inputs.is_empty() {
        return Err(PyValueError::new_err("merge needs at least one input"));
    }
    let inputs_ref: Vec<&str> = inputs.iter().map(|s| s.as_str()).collect();
    merge_native(&inputs_ref, output, force).map_err(|e| PyIOError::new_err(format!("merge: {e}")))
}

/// Pure-Rust entry point used by the rubam-samtools shadow CLI in Phase C.
pub fn merge_native(inputs: &[&str], output: &str, force: bool) -> std::io::Result<()> {
    if !force && std::path::Path::new(output).exists() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::AlreadyExists,
            format!("{output} already exists; pass force=true to overwrite"),
        ));
    }
    if inputs.is_empty() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "merge needs at least one input",
        ));
    }

    // Take the header from the first input.
    let mut first_reader = open_streaming(inputs[0])?;
    let header = read_header_streaming(&mut first_reader)?;

    let mut all_records: Vec<RecordBuf> = Vec::new();
    {
        let mut buf = RecordBuf::default();
        while first_reader.read_record_buf(&header, &mut buf)? != 0 {
            all_records.push(buf.clone());
        }
    }
    for path in &inputs[1..] {
        let mut reader = open_streaming(path)?;
        let _h = read_header_streaming(&mut reader)?;
        let mut buf = RecordBuf::default();
        while reader.read_record_buf(&header, &mut buf)? != 0 {
            all_records.push(buf.clone());
        }
    }
    all_records.sort_by_key(|r| {
        let rid = r.reference_sequence_id().unwrap_or(usize::MAX);
        let pos = r.alignment_start().map(|p| p.get()).unwrap_or(0);
        (rid, pos)
    });

    // Write
    let writer_file = BufWriter::new(File::create(output)?);
    let mut writer: BamWriter = bam::io::Writer::new(writer_file);
    writer.write_header(&header)?;
    for r in &all_records {
        writer.write_alignment_record(&header, r)?;
    }
    writer.try_finish()?;
    Ok(())
}
