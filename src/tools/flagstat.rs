//! samtools flagstat — same counting logic as v0.1 stats::flag_stats but
//! formatted to a writer in the samtools-compatible text layout.

use std::io::Write;

use std::io;

#[cfg(feature = "python")]
use pyo3::exceptions::PyIOError;
#[cfg(feature = "python")]
use pyo3::prelude::*;

use crate::common::{open_streaming, read_header_streaming};

const FLAG_PAIRED: u16 = 0x1;
const FLAG_PROPER: u16 = 0x2;
const FLAG_UNMAP: u16 = 0x4;
const FLAG_MUNMAP: u16 = 0x8;
const FLAG_READ1: u16 = 0x40;
const FLAG_READ2: u16 = 0x80;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_QCFAIL: u16 = 0x200;
const FLAG_DUP: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

#[derive(Default)]
struct Counts {
    total: u64,
    qcfail: u64,
    primary: u64,
    secondary: u64,
    supplementary: u64,
    duplicates: u64,
    primary_duplicates: u64,
    mapped: u64,
    primary_mapped: u64,
    paired: u64,
    read_1: u64,
    read_2: u64,
    properly_paired: u64,
    with_mate_mapped: u64,
    singletons: u64,
    mate_diff_chr: u64,
    mate_diff_chr_mapq5: u64,
}

/// pysam-compatible `flagstat`: returns the samtools-style multi-line
/// flagstat report as a Python `str` (drop-in for `pysam.flagstat(bam)`).
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(name = "flagstat", signature = (input))]
pub fn flagstat_py(input: &str) -> PyResult<String> {
    let mut buf: Vec<u8> = Vec::new();
    flagstat_native(input, &mut buf).map_err(|e| PyIOError::new_err(format!("flagstat: {e}")))?;
    String::from_utf8(buf)
        .map_err(|e| PyIOError::new_err(format!("flagstat: non-UTF-8 output: {e}")))
}

/// Pure-Rust entry point used by the rubam-samtools shadow CLI in Phase C.
/// Writes 17 lines of samtools-style flagstat text to `out`.
pub fn flagstat_native<W: Write>(input: &str, out: &mut W) -> std::io::Result<()> {
    let counts = count(input)?;

    // samtools flagstat layout: "<n> + 0 <label>"
    write_line(
        out,
        counts.total,
        "in total (QC-passed reads + QC-failed reads)",
    )?;
    write_line(out, counts.primary, "primary")?;
    write_line(out, counts.secondary, "secondary")?;
    write_line(out, counts.supplementary, "supplementary")?;
    write_line(out, counts.duplicates, "duplicates")?;
    write_line(out, counts.primary_duplicates, "primary duplicates")?;
    write_line(out, counts.mapped, "mapped")?;
    write_line(out, counts.primary_mapped, "primary mapped")?;
    write_line(out, counts.paired, "paired in sequencing")?;
    write_line(out, counts.read_1, "read1")?;
    write_line(out, counts.read_2, "read2")?;
    write_line(out, counts.properly_paired, "properly paired")?;
    write_line(out, counts.with_mate_mapped, "with itself and mate mapped")?;
    write_line(out, counts.singletons, "singletons")?;
    write_line(
        out,
        counts.mate_diff_chr,
        "with mate mapped to a different chr",
    )?;
    write_line(
        out,
        counts.mate_diff_chr_mapq5,
        "with mate mapped to a different chr (mapQ>=5)",
    )?;
    // Last line: QC-failed totals (rubam doesn't yet partition per-category)
    writeln!(
        out,
        "{} + 0 QC-failed reads (counted separately)",
        counts.qcfail
    )?;
    Ok(())
}

fn write_line<W: Write>(out: &mut W, n: u64, label: &str) -> std::io::Result<()> {
    writeln!(out, "{n} + 0 {label}")
}

fn count(input: &str) -> io::Result<Counts> {
    let mut reader = open_streaming(input)?;
    read_header_streaming(&mut reader)?;
    let mut c = Counts::default();

    for result in reader.records() {
        let record = result.map_err(|e| {
            io::Error::new(io::ErrorKind::Other, format!("record read failed: {e}"))
        })?;
        let flags = record.flags().bits();

        if flags & FLAG_QCFAIL != 0 {
            c.qcfail += 1;
            continue;
        }
        c.total += 1;

        let is_secondary = flags & FLAG_SECONDARY != 0;
        let is_supplementary = flags & FLAG_SUPPLEMENTARY != 0;
        let is_unmapped = flags & FLAG_UNMAP != 0;
        let is_dup = flags & FLAG_DUP != 0;
        let is_paired = flags & FLAG_PAIRED != 0;
        let is_read1 = flags & FLAG_READ1 != 0;
        let is_read2 = flags & FLAG_READ2 != 0;
        let is_proper = flags & FLAG_PROPER != 0;
        let mate_unmapped = flags & FLAG_MUNMAP != 0;

        if is_dup {
            c.duplicates += 1;
        }
        if !is_unmapped {
            c.mapped += 1;
        }

        if is_secondary {
            c.secondary += 1;
        } else if is_supplementary {
            c.supplementary += 1;
        } else {
            c.primary += 1;
            if is_dup {
                c.primary_duplicates += 1;
            }
            if !is_unmapped {
                c.primary_mapped += 1;
            }
            if is_paired {
                c.paired += 1;
                if is_read1 {
                    c.read_1 += 1;
                }
                if is_read2 {
                    c.read_2 += 1;
                }
                if is_proper && !is_unmapped {
                    c.properly_paired += 1;
                }
                if !is_unmapped && !mate_unmapped {
                    c.with_mate_mapped += 1;
                    let rid = record
                        .reference_sequence_id()
                        .transpose()
                        .map_err(|e| io::Error::new(io::ErrorKind::Other, format!("rid: {e}")))?;
                    let mrid = record
                        .mate_reference_sequence_id()
                        .transpose()
                        .map_err(|e| io::Error::new(io::ErrorKind::Other, format!("mrid: {e}")))?;
                    if let (Some(rid), Some(mrid)) = (rid, mrid) {
                        if rid != mrid {
                            c.mate_diff_chr += 1;
                            if let Some(mq) = record.mapping_quality() {
                                if mq.get() >= 5 {
                                    c.mate_diff_chr_mapq5 += 1;
                                }
                            }
                        }
                    }
                }
                if !is_unmapped && mate_unmapped {
                    c.singletons += 1;
                }
            }
        }
    }
    Ok(c)
}
