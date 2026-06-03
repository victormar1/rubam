//! Property-based tests for rubam's CIGAR primitives.
//!
//! These codify the SAMv1 CIGAR consumption rules as proptest invariants so
//! they are continuously checked, not just stated in the manuscript:
//!
//!   - `reference_span(cigar) == sum of M / D / N / = / X op lengths`
//!   - `query_span(cigar)     == sum of M / I / S / = / X op lengths`
//!   - `H` and `P` contribute 0 to both spans
//!   - an empty CIGAR has zero ref-span and zero query-span
//!
//! rubam exposes `api::Cigar::consumes_reference()` / `consumes_query()` on
//! individual ops; the spans are sums of `op.len()` filtered by those
//! predicates. This file proves those predicates correctly identify the SAM
//! consumption sets on arbitrary CIGAR vectors.
//!
//! Generators are bounded (<=50 ops/case, op-len <= 1000) so the full 256-case
//! default suite finishes well under 30 s.

use proptest::prelude::*;
use rubam::api::Cigar;

/// Wall-clock budget guard: keep per-op length small so summed spans stay in
/// `u64` range and the suite stays fast.
const MAX_OP_LEN: u32 = 1_000;
const MAX_OPS: usize = 50;

/// Strategy producing a single arbitrary `Cigar` op, uniform over all 9 kinds
/// with a bounded length.
fn arb_op() -> impl Strategy<Value = Cigar> {
    (0u8..9, 0u32..=MAX_OP_LEN).prop_map(|(tag, n)| match tag {
        0 => Cigar::Match(n),
        1 => Cigar::Ins(n),
        2 => Cigar::Del(n),
        3 => Cigar::RefSkip(n),
        4 => Cigar::Equal(n),
        5 => Cigar::Diff(n),
        6 => Cigar::SoftClip(n),
        7 => Cigar::HardClip(n),
        _ => Cigar::Pad(n),
    })
}

/// Strategy producing an arbitrary CIGAR string of up to `MAX_OPS` ops.
fn arb_cigar() -> impl Strategy<Value = Vec<Cigar>> {
    prop::collection::vec(arb_op(), 0..=MAX_OPS)
}

/// Reference span as defined by SAMv1: sum of lengths of ops that consume the
/// reference (M, D, N, =, X).
fn reference_span(cigar: &[Cigar]) -> u64 {
    cigar
        .iter()
        .filter(|op| op.consumes_reference())
        .map(|op| op.len() as u64)
        .sum()
}

/// Query span as defined by SAMv1: sum of lengths of ops that consume the
/// query (M, I, S, =, X).
fn query_span(cigar: &[Cigar]) -> u64 {
    cigar
        .iter()
        .filter(|op| op.consumes_query())
        .map(|op| op.len() as u64)
        .sum()
}

/// Independent reference computation that pattern-matches each variant
/// directly (no dependency on `consumes_*`). The proptests below check that
/// rubam's `consumes_*` predicates agree with this oracle.
fn reference_span_oracle(cigar: &[Cigar]) -> u64 {
    cigar
        .iter()
        .map(|op| match *op {
            Cigar::Match(n)
            | Cigar::Del(n)
            | Cigar::RefSkip(n)
            | Cigar::Equal(n)
            | Cigar::Diff(n) => n as u64,
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => 0,
        })
        .sum()
}

fn query_span_oracle(cigar: &[Cigar]) -> u64 {
    cigar
        .iter()
        .map(|op| match *op {
            Cigar::Match(n)
            | Cigar::Ins(n)
            | Cigar::SoftClip(n)
            | Cigar::Equal(n)
            | Cigar::Diff(n) => n as u64,
            Cigar::Del(_) | Cigar::RefSkip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => 0,
        })
        .sum()
}

proptest! {
    /// Invariant: rubam's `consumes_reference`-filtered sum agrees with a
    /// hand-written oracle on the M / D / N / = / X consumption set.
    #[test]
    fn reference_span_matches_oracle(cigar in arb_cigar()) {
        prop_assert_eq!(reference_span(&cigar), reference_span_oracle(&cigar));
    }

    /// Invariant: rubam's `consumes_query`-filtered sum agrees with a
    /// hand-written oracle on the M / I / S / = / X consumption set.
    #[test]
    fn query_span_matches_oracle(cigar in arb_cigar()) {
        prop_assert_eq!(query_span(&cigar), query_span_oracle(&cigar));
    }

    /// Invariant: `H` and `P` ops contribute exactly 0 to both spans,
    /// regardless of their length or position in the CIGAR.
    #[test]
    fn hard_clip_and_pad_contribute_zero(
        prefix in arb_cigar(),
        suffix in arb_cigar(),
        h_len in 0u32..=MAX_OP_LEN,
        p_len in 0u32..=MAX_OP_LEN,
    ) {
        let mut without = prefix.clone();
        without.extend(suffix.iter().copied());

        let mut with = prefix;
        with.push(Cigar::HardClip(h_len));
        with.push(Cigar::Pad(p_len));
        with.extend(suffix);

        prop_assert_eq!(reference_span(&without), reference_span(&with));
        prop_assert_eq!(query_span(&without),     query_span(&with));
    }

    /// Invariant: a CIGAR composed of only `H` and `P` ops has zero spans.
    #[test]
    fn only_hard_clip_and_pad_yields_zero_spans(
        ops in prop::collection::vec(
            prop_oneof![
                (0u32..=MAX_OP_LEN).prop_map(Cigar::HardClip),
                (0u32..=MAX_OP_LEN).prop_map(Cigar::Pad),
            ],
            0..=MAX_OPS,
        )
    ) {
        prop_assert_eq!(reference_span(&ops), 0);
        prop_assert_eq!(query_span(&ops),     0);
    }

    /// Invariant: each op is in *exactly* the consumption sets prescribed by
    /// SAMv1. This pins down the predicate tables on every variant.
    #[test]
    fn consumption_predicates_match_samv1(op in arb_op()) {
        let (cr_expected, cq_expected) = match op {
            Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => (true,  true),
            Cigar::Ins(_)   | Cigar::SoftClip(_)               => (false, true),
            Cigar::Del(_)   | Cigar::RefSkip(_)                => (true,  false),
            Cigar::HardClip(_) | Cigar::Pad(_)                 => (false, false),
        };
        prop_assert_eq!(op.consumes_reference(), cr_expected);
        prop_assert_eq!(op.consumes_query(),     cq_expected);
    }
}

#[test]
fn empty_cigar_has_zero_spans() {
    let cigar: Vec<Cigar> = Vec::new();
    assert_eq!(reference_span(&cigar), 0);
    assert_eq!(query_span(&cigar), 0);
    assert_eq!(reference_span_oracle(&cigar), 0);
    assert_eq!(query_span_oracle(&cigar), 0);
}
