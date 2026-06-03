//! Per-base depth: single-region (`get_depths`) and batched (`get_depths_regions`).

use std::sync::Arc;

use noodles::core::{Position, Region};
use noodles::sam::alignment::record::cigar::op::Kind;

#[cfg(feature = "python")]
use numpy::{IntoPyArray, PyArray1};
#[cfg(feature = "python")]
use pyo3::exceptions::{PyIOError, PyValueError};
#[cfg(feature = "python")]
use pyo3::prelude::*;
use rayon::prelude::*;

use crate::common::FLAG_FILTER_DEFAULT;
#[cfg(feature = "python")]
use crate::common::{open_indexed, read_header_indexed, validate_chrom};

/// Memory-vs-throughput tradeoff hint for `get_depths`.
///
/// Default ("fast") preserves v0.3.x behaviour: use every worker the caller
/// requests, no chunk-size cap. The other modes cap parallelism and chunk
/// size to lower peak RSS at a wall-clock cost.
///
/// Rough behaviour on real-WGS chr20 (~64 Mb) measured in Wave 4:
///
/// | mode      | max workers | chunk cap | RSS hint  | speed hint |
/// |-----------|-------------|-----------|-----------|------------|
/// | fast      | num_threads | n/nt      | ≈ 1.0 GB  | 1×         |
/// | balanced  | min(nt, 4)  | 10 Mb     | ≈ 0.8 GB  | ≈ 0.7×     |
/// | low_mem   | min(nt, 2)  | 2 Mb      | ≈ 0.7 GB  | ≈ 0.3×     |
/// | auto      | aliased to "balanced" until v0.5 ships sysinfo detection    |
///
/// Caveat: the dominant cost is the returned `Vec<u64> positions` (~512 MB on
/// chr20). The mode knob reduces per-worker transient buffers, not that final
/// allocation. v0.5 will introduce a `get_depths_only` API that omits the
/// positions vector entirely for callers that only need depths.
#[cfg(feature = "python")]
fn resolve_mode(
    memory_mode: Option<&str>,
    num_threads: usize,
    n: usize,
) -> PyResult<(usize, usize)> // (effective_nt, chunk_size_cap)
{
    let mode = memory_mode.unwrap_or("fast");
    let nt_user = num_threads.max(1).min(n);
    match mode {
        "fast" => Ok((nt_user, n)),  // no cap
        "balanced" | "auto" => Ok((nt_user.min(4), 10_000_000)),
        "low_mem" => Ok((nt_user.min(2),  2_000_000)),
        other => Err(PyValueError::new_err(format!(
            "memory_mode must be one of {{\"fast\", \"balanced\", \"low_mem\", \"auto\"}}, got {other:?}"
        ))),
    }
}

/// Inner depth computation (no pyo3 conversion). Shared between the PyList
/// returning `get_depths` and the numpy returning `get_depths_numpy`.
///
/// **Memory note**: the returned `Vec<u64>` and `Vec<u32>` carry the full
/// per-position result. Converting these to Python via PyList costs ~28 bytes
/// per Python int, dominating peak RSS at ~3.4 GB on chr20 64 Mb. The numpy
/// path materialises them as `np.uint64` / `np.uint32` arrays at 8/4 bytes per
/// element (~768 MB total on chr20), a ~4.5× reduction. Callers that only need
/// depths should pass `step > 1` or use `get_depths_numpy`.
#[cfg(feature = "python")]
#[allow(clippy::too_many_arguments)]
fn compute_depths_internal(
    bam_path: &std::path::Path,
    chromosome: &str,
    start: u64,
    end: u64,
    step: u64,
    min_mapq: u8,
    min_bq: u8,
    max_depth: usize,
    num_threads: usize,
    memory_mode: Option<&str>,
) -> PyResult<(Vec<u64>, Vec<u32>)> {
    let bam_path: &str = bam_path.to_str().ok_or_else(|| {
        PyValueError::new_err(
            "bam_path must be valid UTF-8 (non-UTF-16 Windows paths not yet supported)",
        )
    })?;
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
        return Ok((Vec::new(), Vec::new()));
    }
    let (nt, chunk_cap) = resolve_mode(memory_mode, num_threads, n)?;
    let nt = nt.max(1).min(n);
    // chunk_size bounded by mode cap; n_chunks may exceed nt so the rayon pool
    // serialises excess chunks into fewer in-flight buffers (peak RSS reduction).
    let natural_chunk = n.div_ceil(nt);
    let chunk_size = natural_chunk.min(chunk_cap.max(1));
    let n_chunks = n.div_ceil(chunk_size);

    let header = Arc::new(header);
    let chrom_bytes: Arc<[u8]> = Arc::from(chromosome.as_bytes().to_vec());
    let bam_path: Arc<str> = Arc::from(bam_path);

    let chunk_op = |i: usize| -> Result<Vec<u32>, PyErr> {
        let lo = i * chunk_size;
        let hi = ((i + 1) * chunk_size).min(n);
        if lo >= hi {
            return Ok(Vec::<u32>::new());
        }
        process_depth_chunk(
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
        )
    };

    // Install a bounded rayon pool so peak in-flight chunks <= nt. The default
    // global pool may have more threads available; capping here is what makes
    // memory_mode actually reduce RSS rather than just rename a knob.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nt)
        .build()
        .map_err(|e| PyIOError::new_err(format!("rayon pool: {e}")))?;
    let chunk_results: Result<Vec<Vec<u32>>, PyErr> =
        pool.install(|| (0..n_chunks).into_par_iter().map(chunk_op).collect());

    let depths: Vec<u32> = chunk_results?.into_iter().flatten().collect();
    debug_assert_eq!(depths.len(), positions.len());
    Ok((positions, depths))
}

/// Compute per-base depth over a 1-based, inclusive region. Returns Python
/// `list[int]` for both positions and depths (legacy v0.3.x compatible).
///
/// **For large regions, prefer `get_depths_numpy`** which returns
/// `np.uint64` / `np.uint32` arrays at ~4.5× lower peak RSS than the PyList
/// path on the chr20 64 Mb benchmark.
#[cfg(feature = "python")]
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
    memory_mode = None,
))]
#[allow(clippy::too_many_arguments)]
pub fn get_depths(
    bam_path: std::path::PathBuf,
    chromosome: &str,
    start: u64,
    end: u64,
    step: u64,
    min_mapq: u8,
    min_bq: u8,
    max_depth: usize,
    num_threads: usize,
    memory_mode: Option<&str>,
) -> PyResult<(Vec<u64>, Vec<u32>)> {
    compute_depths_internal(
        bam_path.as_path(),
        chromosome,
        start,
        end,
        step,
        min_mapq,
        min_bq,
        max_depth,
        num_threads,
        memory_mode,
    )
}

/// Numpy-backed variant of `get_depths`. Returns
/// `(np.ndarray[uint64], np.ndarray[uint32])` instead of `(list, list)`.
/// This is the **recommended** entry point for callers that materialise the
/// full positions / depths vectors in-process; PyList materialisation is ~7×
/// more expensive per element (28 bytes per Python int vs 8/4 bytes per numpy
/// element).
///
/// On the HG002 GIAB chr20 64 Mb benchmark (8 threads, fast mode), peak RSS
/// drops from ~3.4 GB (`get_depths`) to ~0.8 GB (`get_depths_numpy`); the
/// numpy arrays are zero-copy views into the underlying Rust `Vec`s and
/// require no element-wise conversion.
#[cfg(feature = "python")]
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
    memory_mode = None,
))]
#[allow(clippy::too_many_arguments)]
pub fn get_depths_numpy<'py>(
    py: Python<'py>,
    bam_path: std::path::PathBuf,
    chromosome: &str,
    start: u64,
    end: u64,
    step: u64,
    min_mapq: u8,
    min_bq: u8,
    max_depth: usize,
    num_threads: usize,
    memory_mode: Option<&str>,
) -> PyResult<(Bound<'py, PyArray1<u64>>, Bound<'py, PyArray1<u32>>)> {
    let (positions, depths) = compute_depths_internal(
        bam_path.as_path(),
        chromosome,
        start,
        end,
        step,
        min_mapq,
        min_bq,
        max_depth,
        num_threads,
        memory_mode,
    )?;
    Ok((positions.into_pyarray(py), depths.into_pyarray(py)))
}

/// Batch version: compute depth for many regions, parallelizing across regions.
///
/// Returns one `(positions, depths)` per region, in the same order as input.
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (
    bam_path,
    regions,
    step = 1,
    min_mapq = 0,
    min_bq = 13,
    max_depth = 8000,
    num_threads = 12,
))]
#[allow(clippy::too_many_arguments)]
pub fn get_depths_regions(
    bam_path: std::path::PathBuf,
    regions: Vec<(String, u64, u64)>,
    step: u64,
    min_mapq: u8,
    min_bq: u8,
    max_depth: usize,
    num_threads: usize,
) -> PyResult<Vec<(Vec<u64>, Vec<u32>)>> {
    if step == 0 {
        return Err(PyValueError::new_err("step must be >= 1"));
    }
    let bam_path: &str = bam_path
        .to_str()
        .ok_or_else(|| PyValueError::new_err("bam_path must be valid UTF-8"))?;

    let mut probe = open_indexed(bam_path)?;
    let header = read_header_indexed(&mut probe)?;
    drop(probe);
    for (chrom, start, end) in &regions {
        if *start == 0 || *end < *start {
            return Err(PyValueError::new_err(format!(
                "invalid region {chrom}:{start}-{end}"
            )));
        }
        validate_chrom(&header, chrom)?;
    }

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads.max(1))
        .build()
        .map_err(|e| PyIOError::new_err(format!("rayon pool: {e}")))?;

    let header = Arc::new(header);
    let bam_path: Arc<str> = Arc::from(bam_path);

    let max_depth_u32 = max_depth as u32;
    let header_ref = header.clone();
    let bam_path_ref = bam_path.clone();

    let results: Result<Vec<(Vec<u64>, Vec<u32>)>, PyErr> = pool.install(|| {
        regions
            .par_iter()
            .map(|(chrom, start, end)| {
                let positions: Vec<u64> = (*start..=*end).step_by(step as usize).collect();
                if positions.is_empty() {
                    return Ok((Vec::<u64>::new(), Vec::<u32>::new()));
                }
                let chrom_bytes: Vec<u8> = chrom.as_bytes().to_vec();
                let depths = process_depth_chunk(
                    &bam_path_ref,
                    &header_ref,
                    &chrom_bytes,
                    positions[0],
                    *positions.last().unwrap(),
                    *start,
                    step,
                    min_mapq,
                    min_bq,
                    max_depth_u32,
                    positions.len(),
                )?;
                Ok((positions, depths))
            })
            .collect()
    });

    results
}

/// Pure-Rust depth computation, no pyo3 dependency. Used by the
/// `rubam-depth` binary for fair benchmarks where every byte of output is
/// produced by Rust (no Python in the hot path).
#[allow(clippy::too_many_arguments)]
pub fn compute_depths_native(
    bam_path: &str,
    chromosome: &str,
    start: u64,
    end: u64,
    step: u64,
    min_mapq: u8,
    min_bq: u8,
    max_depth: u32,
    num_threads: usize,
) -> Result<(Vec<u64>, Vec<u32>), String> {
    if start == 0 || end < start {
        return Err(format!(
            "invalid interval: start={start}, end={end} (require 1<=start<=end)"
        ));
    }
    if step == 0 {
        return Err("step must be >= 1".into());
    }

    let mut probe = noodles::bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)
        .map_err(|e| format!("failed to open indexed BAM at {bam_path}: {e}"))?;
    let header =
        crate::common::read_header_indexed(&mut probe).map_err(|e| format!("read header: {e}"))?;
    if !header
        .reference_sequences()
        .contains_key(chromosome.as_bytes())
    {
        return Err(format!("chromosome {chromosome} not found in BAM header"));
    }
    drop(probe);

    let positions: Vec<u64> = (start..=end).step_by(step as usize).collect();
    let n = positions.len();
    if n == 0 {
        return Ok((Vec::new(), Vec::new()));
    }
    let nt = num_threads.max(1).min(n);
    let chunk_size = n.div_ceil(nt);

    let header = Arc::new(header);
    let chrom_bytes: Arc<[u8]> = Arc::from(chromosome.as_bytes().to_vec());
    let bam_path: Arc<str> = Arc::from(bam_path);

    let chunk_results: Result<Vec<Vec<u32>>, String> = (0..nt)
        .into_par_iter()
        .map(|i| {
            let lo = i * chunk_size;
            let hi = ((i + 1) * chunk_size).min(n);
            if lo >= hi {
                return Ok(Vec::<u32>::new());
            }
            process_depth_chunk_native(
                &bam_path,
                &header,
                &chrom_bytes,
                positions[lo],
                positions[hi - 1],
                start,
                step,
                min_mapq,
                min_bq,
                max_depth,
                hi - lo,
            )
        })
        .collect();

    let depths: Vec<u32> = chunk_results?.into_iter().flatten().collect();
    Ok((positions, depths))
}

#[allow(clippy::too_many_arguments)]
fn process_depth_chunk_native(
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
) -> Result<Vec<u32>, String> {
    let mut reader = noodles::bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)
        .map_err(|e| format!("thread open: {e}"))?;
    crate::common::read_header_indexed(&mut reader)
        .map_err(|e| format!("thread read_header: {e}"))?;

    let region = Region::new(
        chrom.to_vec(),
        Position::new(chunk_start_pos as usize)
            .ok_or_else(|| "chunk_start_pos must be >= 1".to_string())?
            ..=Position::new(chunk_end_pos as usize)
                .ok_or_else(|| "chunk_end_pos must be >= 1".to_string())?,
    );
    let query = reader
        .query(header, &region)
        .map_err(|e| format!("query: {e}"))?;
    let mut depths = vec![0u32; chunk_len];
    for result in query.records() {
        let record = result.map_err(|e| format!("record: {e}"))?;
        if record.flags().bits() & FLAG_FILTER_DEFAULT != 0 {
            continue;
        }
        let mapq_val = record.mapping_quality().map(|q| q.get()).unwrap_or(0);
        if mapq_val < min_mapq {
            continue;
        }
        let aln_start_pos = match record.alignment_start() {
            Some(Ok(p)) => p,
            Some(Err(e)) => return Err(format!("alignment_start: {e}")),
            None => continue,
        };
        let mut ref_pos: u64 = aln_start_pos.get() as u64;
        let mut query_pos: usize = 0;
        let qbytes: Vec<u8> = record.quality_scores().iter().collect();
        let cigar = record.cigar();
        for op_result in cigar.iter() {
            let op = op_result.map_err(|e| format!("cigar: {e}"))?;
            let len = op.len();
            match op.kind() {
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
                        if local_idx < depths.len() && depths[local_idx] < max_depth {
                            depths[local_idx] += 1;
                        }
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
    Ok(depths)
}

#[cfg(feature = "python")]
#[allow(clippy::too_many_arguments)]
fn process_depth_chunk(
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
) -> PyResult<Vec<u32>> {
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

    let mut depths = vec![0u32; chunk_len];

    for result in query.records() {
        let record = result.map_err(|e| PyIOError::new_err(format!("record read failed: {e}")))?;

        if record.flags().bits() & FLAG_FILTER_DEFAULT != 0 {
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
                        if local_idx < depths.len() && depths[local_idx] < max_depth {
                            depths[local_idx] += 1;
                        }
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

    Ok(depths)
}
