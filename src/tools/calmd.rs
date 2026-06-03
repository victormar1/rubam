//! samtools calmd — recompute NM (and, in v0.2.x, MD) by walking each
//! read's CIGAR. v0.2 only updates NM (= total I + total D bases).

#[cfg(feature = "python")]
use std::fs::File;
#[cfg(feature = "python")]
use std::io::BufWriter;

use noodles::bam;
use noodles::bgzf;
use noodles::sam::{
    alignment::io::Write as _, alignment::record::cigar::op::Kind,
    alignment::record::data::field::Tag, alignment::record_buf::data::field::Value,
    alignment::RecordBuf,
};
#[cfg(feature = "python")]
use pyo3::exceptions::PyIOError;
#[cfg(feature = "python")]
use pyo3::prelude::*;

use crate::common::{open_streaming, read_header_streaming};

#[cfg(feature = "python")]
type BamWriter = bam::io::Writer<bgzf::io::Writer<BufWriter<File>>>;

/// `samtools calmd` port. v0.2: writes the input records to `output` with
/// `NM:i:<I+D>` set/replaced. MD reconstruction lands in v0.2.x.
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (input, reference, *, output))]
pub fn calmd(input: &str, reference: &str, output: &str) -> PyResult<()> {
    let f = File::create(output).map_err(|e| PyIOError::new_err(format!("create out: {e}")))?;
    let mut writer: BamWriter = bam::io::Writer::new(BufWriter::new(f));
    calmd_native_to(input, reference, &mut writer)
        .map_err(|e| PyIOError::new_err(format!("calmd: {e}")))
}

/// Pure-Rust entry point used by the rubam-samtools shadow CLI in Phase C.
/// Accepts an open BAM writer so the CLI can route output to stdout.
pub fn calmd_native_to<W: std::io::Write>(
    input: &str,
    reference: &str,
    writer: &mut bam::io::Writer<bgzf::io::Writer<W>>,
) -> std::io::Result<()> {
    let _ = reference; // accepted for samtools API parity; v0.2.x will read it.
    let mut reader = open_streaming(input)?;
    let header = read_header_streaming(&mut reader)?;
    writer.write_header(&header)?;

    let mut buf = RecordBuf::default();
    while reader.read_record_buf(&header, &mut buf)? != 0 {
        // Compute NM = total I + D length.
        let mut nm: u32 = 0;
        for op in buf.cigar().as_ref().iter() {
            match op.kind() {
                Kind::Insertion | Kind::Deletion => nm += op.len() as u32,
                _ => {}
            }
        }
        // Replace or insert NM:i:<nm>. `Data::insert` overwrites on collision
        // (returns the previous (Tag, Value) when present), so this is the
        // correct upsert path. `Value::from(u32)` picks the smallest fitting
        // unsigned variant.
        buf.data_mut().insert(Tag::EDIT_DISTANCE, Value::from(nm));

        writer.write_alignment_record(&header, &buf)?;
    }
    writer.try_finish()?;
    Ok(())
}
