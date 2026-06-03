//! samtools sort — coordinate-sort a BAM by (refid, pos).
//!
//! v0.2 uses an in-memory sort (acceptable for files up to a few hundred MB).
//! External merge-sort with chunked spill files lands in v0.2.x.

use std::fs::File;
use std::io::BufWriter;

use noodles::bam;
use noodles::sam::{
    self, alignment::io::Write as _, alignment::RecordBuf,
    header::record::value::map::header::sort_order::COORDINATE,
    header::record::value::map::header::tag::SORT_ORDER, header::record::value::Map,
};
#[cfg(feature = "python")]
use pyo3::exceptions::PyIOError;
#[cfg(feature = "python")]
use pyo3::prelude::*;

use crate::common::{open_streaming, read_header_streaming};

/// Coordinate-sort a BAM. Equivalent to `samtools sort INPUT -o OUTPUT`.
///
/// `threads` is accepted for API parity with samtools, but v0.2 sorts
/// in a single thread; rayon-driven external merge-sort lands in v0.2.x.
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (input, output, *, threads = 1))]
pub fn sort(input: &str, output: &str, threads: usize) -> PyResult<()> {
    let _ = threads;
    sort_native(input, output).map_err(|e| PyIOError::new_err(format!("sort: {e}")))
}

/// Pure-Rust entry point used by the rubam-samtools shadow CLI in Phase C.
pub fn sort_native(input: &str, output: &str) -> std::io::Result<()> {
    let mut reader = open_streaming(input)?;
    let header = read_header_streaming(&mut reader)?;

    // Materialize records into RecordBuf so we can sort.
    let mut records: Vec<RecordBuf> = Vec::new();
    let mut buf = RecordBuf::default();
    while reader.read_record_buf(&header, &mut buf)? != 0 {
        records.push(buf.clone());
    }
    records.sort_by_key(|r| {
        let rid = r.reference_sequence_id().unwrap_or(usize::MAX);
        let pos = r.alignment_start().map(|p| p.get()).unwrap_or(0);
        (rid, pos)
    });

    // Mark the header SO:coordinate.
    let mut new_header = header.clone();
    *new_header.header_mut() = Some(
        Map::<sam::header::record::value::map::Header>::builder()
            .insert(SORT_ORDER, COORDINATE)
            .build()
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?,
    );

    let writer_file = BufWriter::new(File::create(output)?);
    let mut writer = bam::io::Writer::new(writer_file);
    writer.write_header(&new_header)?;
    for record in &records {
        writer.write_alignment_record(&new_header, record)?;
    }
    writer.finish(&new_header)?;
    Ok(())
}
