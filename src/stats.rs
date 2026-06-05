//! BAM-level statistics: `count_reads` (region) and `flag_stats` (whole file).

use noodles::core::{Position, Region};

use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::common::{
    open_indexed, open_streaming, read_header_indexed, read_header_streaming, validate_chrom,
};

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

/// Count records in a region matching the given filters.
///
/// `flag_required` and `flag_filtered` follow the standard SAM convention:
/// a record is kept iff `(record.flags & flag_required) == flag_required`
/// and `(record.flags & flag_filtered) == 0`.
///
/// **Default semantics differ from `AlignmentFile.count()` -- read this.**
/// This free function defaults to `flag_filtered=0x704`, i.e. it *excludes*
/// unmapped (0x4), secondary (0x100), QC-fail (0x200) and duplicate (0x400)
/// records (supplementary records are still counted). That is the
/// `samtools`-style mask. `AlignmentFile.count()`, by contrast, defaults to
/// pysam's `read_callback='nofilter'` and counts *every* record in the
/// region. The two therefore return different totals on the same region by
/// default -- e.g. 7 here vs 8 from the method when a secondary alignment is
/// present. For pysam-equivalent counting, call `count_reads(..., flag_filtered=0)`;
/// for the samtools-style total from the method, pass `read_callback='all'`.
#[pyfunction]
#[pyo3(signature = (
    bam_path,
    chromosome,
    start,
    end,
    min_mapq = 0,
    flag_required = 0,
    flag_filtered = 0x704,
))]
pub fn count_reads(
    bam_path: &str,
    chromosome: &str,
    start: u64,
    end: u64,
    min_mapq: u8,
    flag_required: u16,
    flag_filtered: u16,
) -> PyResult<u64> {
    if start == 0 || end < start {
        return Err(PyValueError::new_err(format!(
            "invalid interval: start={start}, end={end} (require 1<=start<=end)"
        )));
    }

    let mut reader = open_indexed(bam_path)?;
    let header = read_header_indexed(&mut reader)?;
    validate_chrom(&header, chromosome)?;

    let region = Region::new(
        chromosome.as_bytes().to_vec(),
        Position::new(start as usize).ok_or_else(|| PyValueError::new_err("start must be >= 1"))?
            ..=Position::new(end as usize)
                .ok_or_else(|| PyValueError::new_err("end must be >= 1"))?,
    );

    let query = reader
        .query(&header, &region)
        .map_err(|e| PyIOError::new_err(format!("query failed: {e}")))?;

    let mut count: u64 = 0;
    for result in query.records() {
        let record = result.map_err(|e| PyIOError::new_err(format!("record read failed: {e}")))?;
        let flags = record.flags().bits();
        if (flags & flag_required) != flag_required {
            continue;
        }
        if (flags & flag_filtered) != 0 {
            continue;
        }
        let mapq_val = record.mapping_quality().map(|q| q.get()).unwrap_or(0);
        if mapq_val < min_mapq {
            continue;
        }
        count += 1;
    }
    Ok(count)
}

/// `samtools flagstat` equivalent: returns a `dict` with the per-category counts
/// (QC-passed only — QC-failed reads are reported separately under `qcfail`).
#[pyfunction]
pub fn flag_stats(py: Python<'_>, bam_path: &str) -> PyResult<Py<PyDict>> {
    let mut reader = open_streaming(bam_path)?;
    read_header_streaming(&mut reader)?;

    let mut total: u64 = 0;
    let mut qcfail: u64 = 0;
    let mut secondary: u64 = 0;
    let mut supplementary: u64 = 0;
    let mut duplicates: u64 = 0;
    let mut primary: u64 = 0;
    let mut primary_duplicates: u64 = 0;
    let mut mapped: u64 = 0;
    let mut primary_mapped: u64 = 0;
    let mut paired: u64 = 0;
    let mut read_1: u64 = 0;
    let mut read_2: u64 = 0;
    let mut properly_paired: u64 = 0;
    let mut with_mate_mapped: u64 = 0;
    let mut singletons: u64 = 0;
    let mut mate_diff_chr: u64 = 0;
    let mut mate_diff_chr_mapq5: u64 = 0;

    for result in reader.records() {
        let record = result.map_err(|e| PyIOError::new_err(format!("record read failed: {e}")))?;
        let flags = record.flags().bits();

        if flags & FLAG_QCFAIL != 0 {
            qcfail += 1;
            continue;
        }
        total += 1;

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
            duplicates += 1;
        }
        if !is_unmapped {
            mapped += 1;
        }

        if is_secondary {
            secondary += 1;
        } else if is_supplementary {
            supplementary += 1;
        } else {
            primary += 1;
            if is_dup {
                primary_duplicates += 1;
            }
            if !is_unmapped {
                primary_mapped += 1;
            }
            if is_paired {
                paired += 1;
                if is_read1 {
                    read_1 += 1;
                }
                if is_read2 {
                    read_2 += 1;
                }
                if is_proper && !is_unmapped {
                    properly_paired += 1;
                }
                if !is_unmapped && !mate_unmapped {
                    with_mate_mapped += 1;
                    let rid = record
                        .reference_sequence_id()
                        .transpose()
                        .map_err(|e| PyIOError::new_err(format!("rid: {e}")))?;
                    let mrid = record
                        .mate_reference_sequence_id()
                        .transpose()
                        .map_err(|e| PyIOError::new_err(format!("mrid: {e}")))?;
                    if let (Some(rid), Some(mrid)) = (rid, mrid) {
                        if rid != mrid {
                            mate_diff_chr += 1;
                            if let Some(mq) = record.mapping_quality() {
                                if mq.get() >= 5 {
                                    mate_diff_chr_mapq5 += 1;
                                }
                            }
                        }
                    }
                }
                if !is_unmapped && mate_unmapped {
                    singletons += 1;
                }
            }
        }
    }

    let dict = PyDict::new(py);
    dict.set_item("total", total)?;
    dict.set_item("qcfail", qcfail)?;
    dict.set_item("primary", primary)?;
    dict.set_item("secondary", secondary)?;
    dict.set_item("supplementary", supplementary)?;
    dict.set_item("duplicates", duplicates)?;
    dict.set_item("primary_duplicates", primary_duplicates)?;
    dict.set_item("mapped", mapped)?;
    dict.set_item("primary_mapped", primary_mapped)?;
    dict.set_item("paired", paired)?;
    dict.set_item("read_1", read_1)?;
    dict.set_item("read_2", read_2)?;
    dict.set_item("properly_paired", properly_paired)?;
    dict.set_item("with_itself_and_mate_mapped", with_mate_mapped)?;
    dict.set_item("singletons", singletons)?;
    dict.set_item("mate_mapped_to_different_chr", mate_diff_chr)?;
    dict.set_item("mate_mapped_to_different_chr_mapq_5", mate_diff_chr_mapq5)?;
    Ok(dict.into())
}
