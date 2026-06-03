# 🦀 rubam — fast, native, Windows-ready BAM/VCF toolkit

[![CI](https://github.com/victormar1/rubam/actions/workflows/integration.yaml/badge.svg)](https://github.com/victormar1/rubam/actions/workflows/integration.yaml)
[![Wheels](https://github.com/victormar1/rubam/actions/workflows/wheel-smoke-test.yml/badge.svg)](https://github.com/victormar1/rubam/actions/workflows/wheel-smoke-test.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![Python](https://img.shields.io/badge/python-3.8%20–%203.13-blue)
![Platforms](https://img.shields.io/badge/platforms-Linux%20%7C%20macOS%20%7C%20Windows-informational)

`rubam` is a pure-Rust BAM/VCF analysis library with first-class Python bindings. It provides per-base depth, pileup, flag statistics, read counting, VCF/BCF read+write and indexed query — multi-threaded, with **bit-exact** parity against `samtools` and `pysam`, and **native binaries** for Linux, macOS **and Windows** (no WSL, no MSYS2, no `htslib` system install). The core `AlignmentFile` surface (`fetch` / `count` / `count_coverage` / `pileup` / header access) is a drop-in for `pysam` on BAM, validated base-for-base against `pysam` on real hg38 data. CRAM is **experimental**: `AlignmentFile(path, reference_filename=...)` opens any CRAM and reads its header; record decode is panic-guarded and raises a Python error on codecs `noodles-cram` does not yet support, so it never crashes across the FFI boundary.

> Originally forked from [`rustbam`](https://github.com/shahcompbio/rustbam) (Choi *et al.*). `rubam` is now an independent project: pure-Rust backend ([`noodles`](https://github.com/zaeleus/noodles)), expanded API, full cross-platform CI, and a peer-reviewed validation campaign.

## Why rubam

| Capability | `pysam` | `samtools` CLI | `mosdepth` | **`rubam`** |
|---|:--:|:--:|:--:|:--:|
| Native Windows wheel | ❌ | ❌ | ❌ | ✅ |
| Multi-threaded depth | ❌ (GIL) | ⚠ partial | ❌ | ✅ |
| Python API | ✅ | ❌ | ❌ | ✅ |
| CRAM support | ✅ | ✅ | ✅ | ⚠ (skeleton v0.3.1; full decode v0.4) |
| Pure-Rust (no C dep) | ❌ | n/a | ❌ | ✅ |
| `pip install` "just works" on Windows | ❌ | n/a | ❌ | ✅ |

## Speed

### Synthetic 10 Mb chr20, 30× coverage, 3 reps best-of

| Tool | 1 thread | 8 threads | vs pysam @ 8t |
|---|---:|---:|---:|
| **rubam** | **4.14 s** | **1.51 s** | **6.0×** |
| samtools depth | 8.34 s | 5.79 s | 1.6× |
| pysam | 8.95 s | 9.11 s (GIL) | 1.0× |
| mosdepth | 15.32 s | 13.88 s | 0.66× |

### Real WGS — HG002 GIAB 2x250bp full chr20 (64.4 Mb, ~30× coverage)

Scaling sweep at threads {1, 2, 4, 8, 16}, 3 reps best-of (lower = better):

| Tool | 1t | 2t | 4t | 8t | 16t |
|---|---:|---:|---:|---:|---:|
| **rubam** | 60.4 s | 35.3 s | **21.8 s** | **17.1 s** | **17.1 s** |
| samtools depth | 89.7 s | 43.5 s | 44.1 s | 43.8 s | 45.2 s |
| pysam | 109.7 s | 110.7 s | 111.5 s | 109.4 s | 111.1 s (GIL) |
| mosdepth | 36.7 s | 36.8 s | 36.7 s | 36.3 s | 37.1 s |

rubam scales 3.5× from 1 → 8 threads, then **saturates at 8t** (I/O-bound). samtools scales only 1 → 2 threads. pysam and mosdepth are flat. **At 8 threads, rubam beats every competitor**: 6.4× pysam, 2.6× samtools, 2.1× mosdepth.

### PacBio HiFi long reads — HG002 chr20 1-10 Mb

| Tool | 1 thread | 8 threads | vs pysam @ 8t |
|---|---:|---:|---:|
| **rubam** | **5.3 s** | **1.9 s** | **5.6×** |
| samtools depth | 9.7 s | 5.8 s | 1.8× |
| pysam | 10.6 s | 10.5 s | 1.0× |
| mosdepth | 19.0 s | 19.0 s | 0.55× |

→ rubam handles long-read CIGAR (rich D/I/=) without slowdown.

### RNA-seq spliced reads — synthetic chr20 1-10 Mb, 5% reads with `aM bN cM` CIGAR (intron skip)

| Tool | 1 thread | 8 threads | vs pysam @ 8t |
|---|---:|---:|---:|
| **rubam** | **4.6 s** | **2.3 s** | **4.8×** |
| samtools depth | 8.3 s | 5.8 s | 1.9× |
| pysam | 11.0 s | 10.8 s | 1.0× |

→ rubam correctly skips reference-skip ops (`N`) without crashing; throughput is unchanged vs unspliced data. mosdepth not run on spliced data.

All numbers are best-of-3 wall-clock on the datasets named in each table heading.

## Rust API (for downstream crates)

`rubam` is also a publishable Cargo crate. Add it to your `Cargo.toml`:

```toml
[dependencies]
rubam = "0.3.12"
```

…and use the pure-Rust types directly (no Python, no pyo3):

```rust
use rubam::api::{AlignmentFile, Aux};

fn count_reverse_reads(bam_path: &str) -> rubam::api::Result<usize> {
    let mut bam = AlignmentFile::open(bam_path)?;
    let mut n = 0;
    for r in bam.records() {
        if r?.is_reverse() {
            n += 1;
        }
    }
    Ok(n)
}

fn extract_split_reads(bam_path: &str) -> rubam::api::Result<Vec<String>> {
    let mut bam = AlignmentFile::open(bam_path)?;
    let mut sa_tags = Vec::new();
    for r in bam.records() {
        let r = r?;
        if let Ok(Aux::String(s)) = r.aux(b"SA") {
            sa_tags.push(s.to_owned());
        }
    }
    Ok(sa_tags)
}
```

API surface (v0.2.1, stable):

| Type | Methods |
|---|---|
| `AlignmentFile` | `open(path)`, `header()`, `records()` |
| `Header` | `target_count`, `tid2name(tid)`, `target_len(tid)`, `target_names()` |
| `AlignedSegment` | `qname`, `tid`, `pos`, `mapq`, `seq`, `qual` (raw phred), `seq_len`, 12 flag accessors, `cigar()`, `aux(tag)` |
| `Cigar` | enum with `Match/Ins/Del/RefSkip/Equal/Diff/SoftClip/HardClip/Pad`, each `(u32)` |
| `Aux<'a>` | enum with 18 variants (`Char`, `I8`/`U8`/.../`U32`, `Float`/`Double`, `String`, `HexByteArray`, 8 `Array*`) |

Drop-in replacement for `rust_htslib::bam::Reader::from_path` for codebases that iterate linearly. Indexed query (`fetch`) lands in v0.3.x. The pyo3 wrapper classes (`rubam.AlignmentFile` etc.) coexist with `api::*` and share the same `noodles` backend; v0.2.2 will refactor them to delegate to `api::*` directly.

## Correctness

`rubam` is bit-exact against `samtools depth -a` over **5 × 10⁶ positions** across five datasets, including whole-chromosome chr1:

| Dataset | Positions | rubam vs samtools |
|---|---:|:---:|
| Synthetic chr20 30× WGS | 1 000 000 | **0 mismatches** ✅ |
| Synthetic chr20 spliced (5 % CIGAR `N`) | 1 000 000 | **0 mismatches** ✅ |
| HG002 GIAB 2×250bp chr20 | 1 000 000 | **0 mismatches** ✅ |
| HG002 PacBio HiFi chr20 | 1 000 000 | **0 mismatches** ✅ |
| HG002 GIAB 2×250bp **whole chr1** (249 Mb) | 1 000 000 | **0 mismatches** ✅ |
| **Total** | **5 000 000** | **0 / 5 M ✅** |

VCF-side correctness vs `pysam.VariantFile`: **319 349 / 319 349 = 100.00 %** on the GIAB HG002 truth chr1 (319 k records, 13 MB BGZF).

Cross-tool correctness vs system `bcftools`: **100 %** on `view`, `query`, `sort`.

## Install

```bash
pip install rubam
```

Pre-built wheels are published for Linux, macOS and Windows; a single
`abi3` wheel per OS covers CPython 3.8 → 3.13. No `htslib`, no compiler,
no WSL required — `pip install rubam` just works on Windows.

The NumPy return path (`get_depths_numpy`) needs NumPy at runtime:

```bash
pip install rubam[numpy]
```

## Quick start

```python
import rubam

positions, depths = rubam.get_depths(
    "sample.bam", "chr1", 1_000_000, 1_001_000,
    step=1, min_mapq=20, min_bq=20,
    max_depth=8000, num_threads=12,
)
```

CLI:

```bash
rubam depth sample.bam chr1 1000000 1001000 -n 12 -Q 20 -q 20 > depth.tsv
```

## Features

### Stable (v0.1)
- `get_depths(bam, chr, start, end, ...)` — per-base coverage over a 1-based, inclusive region.
- CLI `rubam depth …`.

### Shipped since v0.1.x
- `count_reads(bam, chr, start, end, ...)` — `pysam.AlignmentFile.count` replacement.
- `flag_stats(bam)` — `samtools flagstat` replacement, returning a Python dict.
- `pileup_bases(bam, chr, start, end, ...)` — A/C/G/T counts per position.
- `get_depths_regions(bam, regions)` — batch BED-style regions with shared thread pool.
- `get_depths_numpy(...)` — zero-copy `np.uint64` / `np.uint32` return path
  (~4.5× lower peak RSS than the list path; needs `pip install rubam[numpy]`).

### Roadmap
- ⚠ CRAM full record decode (v0.4): `rubam.AlignmentFile("sample.cram", reference_filename="ref.fa")` already opens and reads the header; record decode is panic-guarded and raises a Python error on codecs `noodles-cram` does not yet support (e.g. Huffman byte-series on NYGC-style CRAMs). Tracking the upstream codec landing.
- `to_pandas()` zero-copy helper; Parquet output.
- `rubam.compat.pysam` drop-in shim (v0.5).

### What's new in 0.2

- `rubam.AlignmentFile` and `rubam.AlignedSegment` — drop-in pysam-style
  read iteration and per-read property access (flags, cigar, sequence,
  qualities, tags, reference helpers).
- `AlignmentFile.fetch(chr, start, end)` — indexed region iterator.
- `AlignmentFile.pileup(chr, start, end)` — buffered per-position iterator
  yielding `PileupColumn` objects with `(reference_pos, depth, A/C/G/T/N)`.
- `rubam.tools.{sort, index, view, merge, flagstat, idxstats, calmd, faidx}`
  — pure-Rust ports of the eight most-used samtools subcommands.
- `rubam-samtools` shadow CLI binary — `alias samtools='rubam samtools'`
  and your shell pipelines keep working, on Windows included.

### What's new in 0.2.1

- `rubam::api::{AlignmentFile, AlignedSegment, Header, Cigar, Aux, Error}`
  — pure-Rust public crate API. External Rust crates drop in
  `rubam = "0.3.12"` and import these types directly without pulling
  in pyo3 — a drop-in for `rust-htslib::bam::Reader` for codebases
  that iterate linearly. The public surface is pinned by
  `tests/api_smoke.rs` and `tests/integration_test.rs`.

### What's new in 0.3

- `rubam.VariantFile` and `rubam.VariantRecord` — pysam-style VCF / BCF /
  Tabix support. Read, write (modes `"w"` / `"wz"` / `"wb"` for plain /
  BGZF / BCF), iterate, indexed `fetch(contig, start, end)`, multi-sample
  genotype access via `record.samples["NA12878"]["GT"]`.
- `rubam.VariantHeader` — read-only metadata: samples, contigs (with
  lengths), INFO / FORMAT meta lines (id / number / type / description),
  FILTER ids, file format version.
- `rubam.VariantRecord(header=, …)` constructor — build records from
  scratch. Plus `set_position`, `set_quality`, `set_filter`, `add_filter`,
  `clear_filters`, `set_info` mutation APIs.
- `rubam.tools.bcftools.{view, norm, concat, query, index, sort, stats}`
  — pure-Rust ports of seven most-used bcftools subcommands.
- `rubam-bcftools` shadow CLI — `alias bcftools='rubam bcftools'` works
  on Windows. Same shape as `rubam-samtools`.
- Cross-tool correctness: `(chrom, pos, ref, alt, ids, qual, filters)`
  extracted via both `rubam.VariantFile` and `pysam.VariantFile` agree on
  **0 / 100 records mismatch on a 3-sample synthetic VCF**.

### What's new in 0.3.12

pysam parity on real-world hg38 BAMs, verified base-for-base against
pysam 0.24.0 (`tests/test_pysam_parity_findings.py`):

- **Tolerant header parsing** — opens real hg38 / GATK / Picard BAMs that
  a strict SAM-header parser rejects (`@HD` with no `VN`, multi-part
  versions like `VN:1.6.0`, duplicate `@PG`/`@RG`/`@SQ` IDs from re-run
  pipelines). Valid headers still take the strict fast path unchanged.
- **`count` matches pysam defaults** — `read_callback='nofilter'` by
  default (counts every read in the region, including secondary /
  supplementary / duplicate / QC-fail); `read_callback='all'` applies the
  `0x704` mask.
- **`count_coverage` matches pysam defaults** — `quality_threshold=15`
  (base counted iff `qual >= threshold`), no depth cap, and a
  `read_callback` argument.

The compatibility layer `rubam.compat.pysam` (drop-in `from rubam.compat import pysam`)
lands in v0.5; v0.2 + v0.3 are the foundation it sits on top of.

## Validation & benchmarks

`rubam` is validated against `pysam`, `samtools depth`, `samtools mpileup`, `mosdepth`, `bedtools genomecov` and the original `rustbam` on real WGS, RNA-seq, exome and PacBio HiFi datasets (HG002, NA12878, public ENA RNA-seq), with multi-threaded scaling and cross-platform parity. The numbers in the tables above are drawn from that campaign.

## License

MIT — see `LICENSE`.

## Citation

If you use `rubam` in academic work, please cite the bioRxiv preprint (link will be added once posted).
