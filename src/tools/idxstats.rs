//! samtools idxstats — per-contig (length, mapped, unmapped) using the BAI metadata.

use std::io::Write;

#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::types::{PyDict, PyList};

use crate::common::{open_indexed, read_header_indexed};

/// `samtools idxstats` equivalent. Returns a list of dicts:
/// `{"contig": str, "length": int, "mapped": int, "unmapped": int}`,
/// one entry per reference sequence in the order they appear in the header.
#[cfg(feature = "python")]
#[pyfunction]
pub fn idxstats<'py>(py: Python<'py>, input: &str) -> PyResult<Bound<'py, PyList>> {
    let mut reader = open_indexed(input)?;
    let header = read_header_indexed(&mut reader)?;
    let index = reader.index();

    let pylist = PyList::empty(py);
    let ref_seqs_iter = index.reference_sequences();
    for ((name, ref_seq), idx_ref) in header.reference_sequences().iter().zip(ref_seqs_iter) {
        let metadata = idx_ref.metadata();
        let mapped = metadata.map(|m| m.mapped_record_count()).unwrap_or(0);
        let unmapped = metadata.map(|m| m.unmapped_record_count()).unwrap_or(0);
        let row = PyDict::new(py);
        row.set_item("contig", String::from_utf8_lossy(name).into_owned())?;
        row.set_item("length", ref_seq.length().get())?;
        row.set_item("mapped", mapped)?;
        row.set_item("unmapped", unmapped)?;
        pylist.append(row)?;
    }
    Ok(pylist)
}

/// Pure-Rust entry point used by the rubam-samtools shadow CLI in Phase C.
/// Writes one tab-separated line per contig to `out`:
/// `chrom\tlength\tmapped\tunmapped\n`.
pub fn idxstats_native<W: Write>(input: &str, out: &mut W) -> std::io::Result<()> {
    let mut reader = open_indexed(input)?;
    let header = read_header_indexed(&mut reader)?;
    let index = reader.index();
    let ref_seqs_iter = index.reference_sequences();
    for ((name, ref_seq), idx_ref) in header.reference_sequences().iter().zip(ref_seqs_iter) {
        let metadata = idx_ref.metadata();
        let mapped = metadata.map(|m| m.mapped_record_count()).unwrap_or(0);
        let unmapped = metadata.map(|m| m.unmapped_record_count()).unwrap_or(0);
        writeln!(
            out,
            "{}\t{}\t{}\t{}",
            String::from_utf8_lossy(name),
            ref_seq.length().get(),
            mapped,
            unmapped,
        )?;
    }
    Ok(())
}
