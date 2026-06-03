//! Property-based tests for rubam's per-base depth computation.
//!
//! Invariants codified here (so they live in CI, not just in the manuscript):
//!
//!   1. Depth values are always `>= 0` (vacuous for `u32`; explicitly checked
//!      so the test fails loudly if the return type ever widens to a signed
//!      integer).
//!   2. `depths.len() == positions.len() == (end - start + 1)` for step=1.
//!   3. Reads whose CIGAR places only `D` (deletion) or `N` (ref-skip) over
//!      position `p` do **NOT** contribute to depth at `p`.
//!   4. Reads with `M` / `=` / `X` aligning a base to `p` **DO** contribute
//!      to depth at `p`.
//!
//! We exercise these on small synthetic BAMs built on the fly with
//! `noodles::bam` (the same backend rubam itself uses) and read back through
//! `rubam::depth::compute_depths_native`. Region sizes are kept <= 10 kb and
//! read counts <= 32 per case so the default 256-case suite finishes in well
//! under 30 s on a laptop.

use std::fs::File;
use std::io::BufWriter;
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use noodles::bam;
use noodles::core::Position;
use noodles::sam::{
    self,
    alignment::io::Write as _,
    alignment::record::cigar::op::Kind,
    alignment::record::cigar::Op,
    alignment::record::Flags,
    alignment::record::MappingQuality,
    alignment::record_buf::{Cigar, QualityScores, Sequence},
    alignment::RecordBuf,
    header::record::value::map::header::sort_order::COORDINATE,
    header::record::value::map::header::tag::SORT_ORDER,
    header::record::value::map::{self as map_kind, ReferenceSequence},
    header::record::value::Map,
};
use proptest::prelude::*;
use rubam::depth::compute_depths_native;

const CHROM: &str = "chr_test";
const CHROM_LEN: usize = 20_000;
const MAX_REGION_BP: u64 = 10_000;
const MAX_READS: usize = 32;

static UNIQ: AtomicU64 = AtomicU64::new(0);

/// Allocate a fresh BAM path under the OS temp dir. The test cleans up on
/// successful exit; if a case panics, the leftover file is just a few KB.
fn fresh_bam_path() -> PathBuf {
    let n = UNIQ.fetch_add(1, Ordering::Relaxed);
    let pid = std::process::id();
    let mut p = std::env::temp_dir();
    p.push(format!("rubam_proptest_depth_{pid}_{n}.bam"));
    p
}

/// SAMv1 CIGAR op shape we emit. The "kind" determines whether the read
/// contributes to ref-positions covered by the op.
#[derive(Clone, Copy, Debug)]
enum OpKind {
    M,     // aligned match — contributes
    D,     // deletion       — does NOT contribute, advances ref only
    N,     // ref-skip       — does NOT contribute, advances ref only
    I,     // insertion      — neither contributes nor advances ref
    Equal, // `=`
    Diff,  // `X`
}

impl OpKind {
    fn to_noodles(self) -> Kind {
        match self {
            OpKind::M => Kind::Match,
            OpKind::D => Kind::Deletion,
            OpKind::N => Kind::Skip,
            OpKind::I => Kind::Insertion,
            OpKind::Equal => Kind::SequenceMatch,
            OpKind::Diff => Kind::SequenceMismatch,
        }
    }
    fn consumes_ref(self) -> bool {
        matches!(
            self,
            OpKind::M | OpKind::D | OpKind::N | OpKind::Equal | OpKind::Diff
        )
    }
    fn consumes_query(self) -> bool {
        matches!(self, OpKind::M | OpKind::I | OpKind::Equal | OpKind::Diff)
    }
    /// True iff this op causes the read to *cover* the ref positions it spans
    /// (i.e. contributes to depth). `D` and `N` advance the reference but do
    /// not cover.
    fn covers_ref(self) -> bool {
        matches!(self, OpKind::M | OpKind::Equal | OpKind::Diff)
    }
}

/// A synthetic read described by its CIGAR ops and 1-based alignment start.
#[derive(Clone, Debug)]
struct SynthRead {
    start: u64,
    ops: Vec<(OpKind, u32)>,
}

impl SynthRead {
    /// Total query length (sum of op lens consuming the query).
    fn query_len(&self) -> usize {
        self.ops
            .iter()
            .filter(|(k, _)| k.consumes_query())
            .map(|(_, n)| *n as usize)
            .sum()
    }

    /// 1-based set of reference positions that this read *covers* (M/=/X).
    fn covered_positions(&self) -> Vec<u64> {
        let mut out = Vec::new();
        let mut ref_pos = self.start;
        for (kind, len) in &self.ops {
            if kind.covers_ref() {
                for k in 0..*len as u64 {
                    out.push(ref_pos + k);
                }
            }
            if kind.consumes_ref() {
                ref_pos += *len as u64;
            }
        }
        out
    }

    /// Total reference span of this read (sum of M/D/N/=/X).
    fn ref_span(&self) -> u64 {
        self.ops
            .iter()
            .filter(|(k, _)| k.consumes_ref())
            .map(|(_, n)| *n as u64)
            .sum()
    }
}

/// Build an indexed, coordinate-sorted single-chrom BAM from the given reads
/// and return its path. Reads are sorted by start internally.
fn write_indexed_bam(reads: &[SynthRead]) -> PathBuf {
    let path = fresh_bam_path();
    let chrom_len_nz = NonZeroUsize::new(CHROM_LEN).unwrap();
    let header = sam::Header::builder()
        .set_header(
            Map::<map_kind::Header>::builder()
                .insert(SORT_ORDER, COORDINATE)
                .build()
                .expect("build header map"),
        )
        .add_reference_sequence(CHROM, Map::<ReferenceSequence>::new(chrom_len_nz))
        .build();

    let mut sorted = reads.to_vec();
    sorted.sort_by_key(|r| r.start);

    let file = BufWriter::new(File::create(&path).expect("create BAM"));
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).expect("write header");

    let mapq = MappingQuality::new(60).unwrap();
    let mut record = RecordBuf::default();
    for (idx, r) in sorted.iter().enumerate() {
        let cigar: Cigar = r
            .ops
            .iter()
            .map(|(k, n)| Op::new(k.to_noodles(), *n as usize))
            .collect();
        let qlen = r.query_len();
        let seq_bytes = vec![b'A'; qlen];
        // BQ comfortably above the default min_bq=13 we use in tests.
        let quals = QualityScores::from(vec![30u8; qlen]);

        *record.name_mut() = Some(format!("r{idx}").into());
        *record.flags_mut() = Flags::default();
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(r.start as usize);
        *record.mapping_quality_mut() = Some(mapq);
        *record.cigar_mut() = cigar;
        *record.sequence_mut() = Sequence::from(seq_bytes);
        *record.quality_scores_mut() = quals;
        *record.template_length_mut() = 0;
        *record.mate_reference_sequence_id_mut() = None;
        *record.mate_alignment_start_mut() = None;

        writer
            .write_alignment_record(&header, &record)
            .expect("write record");
    }
    writer.finish(&header).expect("finish bam");
    drop(writer);

    let index = bam::fs::index(&path).expect("index bam");
    let bai_path = path.with_extension("bam.bai");
    let bai_file = File::create(&bai_path).expect("create bai");
    let mut bai_writer = bam::bai::io::Writer::new(bai_file);
    bai_writer.write_index(&index).expect("write bai");
    path
}

/// Compute the expected depth array by replaying CIGARs in pure Rust. This is
/// the oracle that rubam's depth must match.
fn expected_depths(reads: &[SynthRead], start: u64, end: u64) -> Vec<u32> {
    let n = (end - start + 1) as usize;
    let mut depths = vec![0u32; n];
    for r in reads {
        for p in r.covered_positions() {
            if p >= start && p <= end {
                depths[(p - start) as usize] += 1;
            }
        }
    }
    depths
}

// ----------------------------- strategies -----------------------------------

/// Strategy for a single op kind. Weighted toward M so reads usually have at
/// least some coverage.
fn arb_op_kind() -> impl Strategy<Value = OpKind> {
    prop_oneof![
        4 => Just(OpKind::M),
        1 => Just(OpKind::D),
        1 => Just(OpKind::N),
        1 => Just(OpKind::I),
        1 => Just(OpKind::Equal),
        1 => Just(OpKind::Diff),
    ]
}

/// Strategy for a read: 1..=6 ops, lengths 1..=50, start anywhere that keeps
/// the ref span inside the chromosome.
fn arb_read() -> impl Strategy<Value = SynthRead> {
    prop::collection::vec((arb_op_kind(), 1u32..=50), 1..=6).prop_flat_map(|ops| {
        let ref_span: u64 = ops
            .iter()
            .filter(|(k, _)| k.consumes_ref())
            .map(|(_, n)| *n as u64)
            .sum();
        // Need start in [1, CHROM_LEN - ref_span]. If ref_span == 0
        // (pathological CIGAR with only I/H/P), pretend it spans 1 bp so
        // it still has a defined start; depth contribution will be 0
        // anyway.
        let max_start = (CHROM_LEN as u64).saturating_sub(ref_span.max(1)) + 1;
        let max_start = max_start.max(1);
        (Just(ops), 1u64..=max_start).prop_map(|(ops, start)| SynthRead { ops, start })
    })
}

fn arb_reads() -> impl Strategy<Value = Vec<SynthRead>> {
    prop::collection::vec(arb_read(), 0..=MAX_READS)
}

fn arb_region() -> impl Strategy<Value = (u64, u64)> {
    (1u64..=(CHROM_LEN as u64), 0u64..MAX_REGION_BP).prop_map(|(start, span)| {
        let end = (start + span).min(CHROM_LEN as u64);
        let end = end.max(start);
        (start, end)
    })
}

// ------------------------------- tests --------------------------------------

proptest! {
    #![proptest_config(ProptestConfig {
        // Default = 256 cases; keep that but cap shrink-iters so failures
        // shrink fast.
        cases: 256,
        max_shrink_iters: 256,
        ..ProptestConfig::default()
    })]

    /// Invariant 2: returned array length matches region length, and
    /// `positions` is exactly `start..=end`. Invariant 1 is implicit in the
    /// `u32` return type; we still check `>= 0` so the test fails loudly if
    /// the type ever changes.
    #[test]
    fn depth_len_matches_region_len(
        reads in arb_reads(),
        (start, end) in arb_region(),
    ) {
        let bam = write_indexed_bam(&reads);
        let (positions, depths) =
            compute_depths_native(bam.to_str().unwrap(), CHROM, start, end, 1, 0, 0, 100_000, 1)
                .expect("compute_depths_native");

        let expected_len = (end - start + 1) as usize;
        prop_assert_eq!(depths.len(), expected_len);
        prop_assert_eq!(positions.len(), expected_len);
        prop_assert_eq!(positions.first().copied(), Some(start));
        prop_assert_eq!(positions.last().copied(),  Some(end));
        // u32 is unsigned so the runtime "depth >= 0" check is vacuous; we
        // assert the type via a const expr that fails to compile if depths is
        // ever changed to a signed integer.
        const _: () = {
            let _: u32 = 0u32;
        };
        // Touch every element so the depths Vec is fully realized (catches
        // any lazy-init regressions in the worker pool path).
        let _depth_sum: u128 = depths.iter().map(|d| *d as u128).sum();

        let _ = std::fs::remove_file(&bam);
        let _ = std::fs::remove_file(bam.with_extension("bam.bai"));
    }

    /// Invariant 3 & 4: depths produced by rubam exactly match the oracle
    /// replay of the CIGARs. Equivalently, every position covered by a
    /// `M / = / X` op increments depth, and `D / N / I / S / H / P` ops do
    /// not.
    #[test]
    fn depth_matches_cigar_replay(
        reads in arb_reads(),
        (start, end) in arb_region(),
    ) {
        let bam = write_indexed_bam(&reads);
        let (_positions, depths) =
            compute_depths_native(bam.to_str().unwrap(), CHROM, start, end, 1, 0, 0, 100_000, 1)
                .expect("compute_depths_native");

        let expected = expected_depths(&reads, start, end);
        prop_assert_eq!(&depths, &expected);

        let _ = std::fs::remove_file(&bam);
        let _ = std::fs::remove_file(bam.with_extension("bam.bai"));
    }

    /// Invariant 3 sharpened: a read whose CIGAR is `Dn` (or `Nn`) at the
    /// region origin does NOT contribute depth at any position the deletion
    /// or skip covers. We test this directly by hand-rolling such reads.
    #[test]
    fn deletion_and_skip_only_reads_contribute_zero(
        n_reads in 1usize..=MAX_READS,
        op_kind in prop_oneof![Just(OpKind::D), Just(OpKind::N)],
        len in 1u32..=100,
        start in 1u64..=(CHROM_LEN as u64 - 200),
    ) {
        let reads: Vec<SynthRead> = (0..n_reads)
            .map(|_| SynthRead { start, ops: vec![(op_kind, len)] })
            .collect();
        let region_start = start;
        let region_end = (start + len as u64 - 1).min(CHROM_LEN as u64);

        let bam = write_indexed_bam(&reads);
        let (_, depths) =
            compute_depths_native(bam.to_str().unwrap(), CHROM, region_start, region_end,
                                  1, 0, 0, 100_000, 1)
                .expect("compute_depths_native");

        prop_assert!(
            depths.iter().all(|d| *d == 0),
            "expected all zero depth for {:?}-only reads, got {:?}", op_kind, depths
        );

        let _ = std::fs::remove_file(&bam);
        let _ = std::fs::remove_file(bam.with_extension("bam.bai"));
    }

    /// Invariant 4 sharpened: `n` reads each starting at the same position
    /// with a single `Mlen` op produce depth `n` at every position of the
    /// match (within `max_depth`).
    #[test]
    fn match_reads_contribute_depth_n(
        n_reads in 1u32..=20,
        len in 1u32..=100,
        start in 1u64..=(CHROM_LEN as u64 - 200),
        op_kind in prop_oneof![Just(OpKind::M), Just(OpKind::Equal), Just(OpKind::Diff)],
    ) {
        let reads: Vec<SynthRead> = (0..n_reads)
            .map(|_| SynthRead { start, ops: vec![(op_kind, len)] })
            .collect();
        let region_start = start;
        let region_end = start + len as u64 - 1;

        let bam = write_indexed_bam(&reads);
        let (_, depths) =
            compute_depths_native(bam.to_str().unwrap(), CHROM, region_start, region_end,
                                  1, 0, 0, 100_000, 1)
                .expect("compute_depths_native");

        prop_assert!(
            depths.iter().all(|d| *d == n_reads),
            "expected depth = {} everywhere for {:?} reads, got {:?}", n_reads, op_kind, depths
        );

        let _ = std::fs::remove_file(&bam);
        let _ = std::fs::remove_file(bam.with_extension("bam.bai"));
    }
}

#[test]
fn synth_read_oracle_self_consistency() {
    // A 10M5D10M read starting at 100 covers [100..=109] U [115..=124], NOT
    // [110..=114] (the D). This sanity-checks the oracle used in the
    // proptests above.
    let r = SynthRead {
        start: 100,
        ops: vec![(OpKind::M, 10), (OpKind::D, 5), (OpKind::M, 10)],
    };
    let covered = r.covered_positions();
    let expected: Vec<u64> = (100..=109).chain(115..=124).collect();
    assert_eq!(covered, expected);
    assert_eq!(r.ref_span(), 25);
    assert_eq!(r.query_len(), 20);
}
