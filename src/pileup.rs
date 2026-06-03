//! Per-base A/C/G/T/N counts (`pileup_bases`) — `samtools mpileup`-style coverage by base.

use std::sync::Arc;

use noodles::core::{Position, Region};
use noodles::sam::alignment::record::cigar::op::Kind;

use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use rayon::prelude::*;

use crate::common::{open_indexed, read_header_indexed, validate_chrom};

/// Per-position counts for each nucleotide.
///
/// Returns a 7-tuple `(positions, a, c, g, t, n, depth)` where `depth = a+c+g+t+n`.
/// Each list has the same length as `positions`. `n` includes both `N` bases and
/// any base that does not match `A/C/G/T/N` after upper-casing.
#[pyfunction]
#[pyo3(signature = (
    bam_path,
    chromosome,
    start,
    end,
    step = 1,
    min_mapq = 0,
    min_bq = 13,
    max_depth = 8000,
    num_threads = 12,
    flag_filter = 0x704,
))]
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn pileup_bases(
    bam_path: &str,
    chromosome: &str,
    start: u64,
    end: u64,
    step: u64,
    min_mapq: u8,
    min_bq: u8,
    max_depth: usize,
    num_threads: usize,
    flag_filter: u16,
) -> PyResult<(
    Vec<u64>,
    Vec<u32>,
    Vec<u32>,
    Vec<u32>,
    Vec<u32>,
    Vec<u32>,
    Vec<u32>,
)> {
    if start == 0 || end < start {
        return Err(PyValueError::new_err(format!(
            "invalid interval: start={start}, end={end} (require 1<=start<=end)"
        )));
    }
    if step == 0 {
        return Err(PyValueError::new_err("step must be >= 1"));
    }

    let mut probe = open_indexed(bam_path)?;
    let header = read_header_indexed(&mut probe)?;
    validate_chrom(&header, chromosome)?;
    drop(probe);

    let positions: Vec<u64> = (start..=end).step_by(step as usize).collect();
    let n = positions.len();
    if n == 0 {
        return Ok((
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        ));
    }
    let nt = num_threads.max(1).min(n);
    let chunk_size = n.div_ceil(nt);

    let header = Arc::new(header);
    let chrom_bytes: Arc<[u8]> = Arc::from(chromosome.as_bytes().to_vec());
    let bam_path: Arc<str> = Arc::from(bam_path);

    let chunk_results: Result<Vec<PileupChunk>, PyErr> = (0..nt)
        .into_par_iter()
        .map(|i| {
            let lo = i * chunk_size;
            let hi = ((i + 1) * chunk_size).min(n);
            if lo >= hi {
                return Ok(PileupChunk::default());
            }
            process_pileup_chunk(
                &bam_path,
                &header,
                &chrom_bytes,
                positions[lo],
                positions[hi - 1],
                start,
                step,
                min_mapq,
                min_bq,
                max_depth as u32,
                hi - lo,
                flag_filter,
            )
        })
        .collect();

    let mut a = Vec::with_capacity(n);
    let mut c = Vec::with_capacity(n);
    let mut g = Vec::with_capacity(n);
    let mut t = Vec::with_capacity(n);
    let mut nn = Vec::with_capacity(n);
    let mut depth = Vec::with_capacity(n);
    for chunk in chunk_results? {
        a.extend(chunk.a);
        c.extend(chunk.c);
        g.extend(chunk.g);
        t.extend(chunk.t);
        nn.extend(chunk.n);
        depth.extend(chunk.depth);
    }

    Ok((positions, a, c, g, t, nn, depth))
}

#[derive(Default)]
struct PileupChunk {
    a: Vec<u32>,
    c: Vec<u32>,
    g: Vec<u32>,
    t: Vec<u32>,
    n: Vec<u32>,
    depth: Vec<u32>,
}

#[allow(clippy::too_many_arguments)]
fn process_pileup_chunk(
    bam_path: &str,
    header: &noodles::sam::Header,
    chrom: &[u8],
    chunk_start_pos: u64,
    chunk_end_pos: u64,
    sample_origin: u64,
    step: u64,
    min_mapq: u8,
    min_bq: u8,
    max_depth: u32,
    chunk_len: usize,
    flag_filter: u16,
) -> PyResult<PileupChunk> {
    let mut reader = open_indexed(bam_path)?;
    read_header_indexed(&mut reader)?;

    let region = Region::new(
        chrom.to_vec(),
        Position::new(chunk_start_pos as usize)
            .ok_or_else(|| PyValueError::new_err("chunk_start_pos must be >= 1"))?
            ..=Position::new(chunk_end_pos as usize)
                .ok_or_else(|| PyValueError::new_err("chunk_end_pos must be >= 1"))?,
    );

    let query = reader
        .query(header, &region)
        .map_err(|e| PyIOError::new_err(format!("query failed: {e}")))?;

    let mut a = vec![0u32; chunk_len];
    let mut c = vec![0u32; chunk_len];
    let mut g = vec![0u32; chunk_len];
    let mut t = vec![0u32; chunk_len];
    let mut n_arr = vec![0u32; chunk_len];
    let mut depth = vec![0u32; chunk_len];

    for result in query.records() {
        let record = result.map_err(|e| PyIOError::new_err(format!("record read failed: {e}")))?;

        if record.flags().bits() & flag_filter != 0 {
            continue;
        }
        let mapq_val = record.mapping_quality().map(|q| q.get()).unwrap_or(0);
        if mapq_val < min_mapq {
            continue;
        }

        let aln_start_pos = match record.alignment_start() {
            Some(Ok(p)) => p,
            Some(Err(e)) => {
                return Err(PyIOError::new_err(format!(
                    "alignment_start parse error: {e}"
                )));
            }
            None => continue,
        };
        let mut ref_pos: u64 = aln_start_pos.get() as u64;
        let mut query_pos: usize = 0;

        let qbytes: Vec<u8> = record.quality_scores().iter().collect();
        let sequence = record.sequence();
        let cigar = record.cigar();

        for op_result in cigar.iter() {
            let op =
                op_result.map_err(|e| PyIOError::new_err(format!("cigar parse error: {e}")))?;
            let len = op.len();
            let kind = op.kind();

            match kind {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    let op_ref_end = ref_pos + len as u64;
                    if ref_pos > chunk_end_pos {
                        let _ = op_ref_end;
                        break;
                    }
                    for k in 0..len {
                        let p = ref_pos + k as u64;
                        if p < chunk_start_pos {
                            continue;
                        }
                        if p > chunk_end_pos {
                            break;
                        }
                        if (p - sample_origin) % step != 0 {
                            continue;
                        }
                        let bq = qbytes.get(query_pos + k).copied().unwrap_or(0);
                        if bq < min_bq {
                            continue;
                        }
                        let local_idx = ((p - chunk_start_pos) / step) as usize;
                        if local_idx >= depth.len() || depth[local_idx] >= max_depth {
                            continue;
                        }
                        let base = sequence
                            .get(query_pos + k)
                            .map(|b| b.to_ascii_uppercase())
                            .unwrap_or(b'N');
                        match base {
                            b'A' => a[local_idx] += 1,
                            b'C' => c[local_idx] += 1,
                            b'G' => g[local_idx] += 1,
                            b'T' => t[local_idx] += 1,
                            _ => n_arr[local_idx] += 1,
                        }
                        depth[local_idx] += 1;
                    }
                    ref_pos = op_ref_end;
                    query_pos += len;
                }
                Kind::Deletion | Kind::Skip => {
                    ref_pos += len as u64;
                }
                Kind::Insertion | Kind::SoftClip => {
                    query_pos += len;
                }
                Kind::HardClip | Kind::Pad => {}
            }
        }
    }

    Ok(PileupChunk {
        a,
        c,
        g,
        t,
        n: n_arr,
        depth,
    })
}
