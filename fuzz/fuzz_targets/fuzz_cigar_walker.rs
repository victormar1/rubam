#![no_main]
//! Fuzz target: arbitrary CIGAR strings → noodles parse → walk and compute
//! reference-span and query-span. Mirrors what
//! `rubam::api::cigar::Cigar::consumes_reference` /
//! `Cigar::consumes_query` are used for in callers like HARMOS, even though
//! the helper functions themselves don't live in `src/api/cigar.rs` (they
//! live inline in the depth/pileup walkers — see `src/depth.rs` and
//! `src/pileup.rs`). The fuzzer here re-implements the walk against the
//! same op semantics so a panic in the noodles CIGAR parser or in our
//! length arithmetic is caught.

use libfuzzer_sys::fuzz_target;

use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::record::Cigar as TextCigar;

fuzz_target!(|data: &[u8]| {
    // Bound the input length to keep iterations fast.
    if data.len() > 4096 {
        return;
    }
    // Feed raw bytes to the noodles SAM-text CIGAR parser. Most inputs
    // get rejected by `parse_op` (per-iteration Result::Err), never a panic.
    let cigar = TextCigar::new(data);

    // Walk the ops. Each `op` is a Result<Op, _> in noodles' iterator.
    let mut ref_span: u64 = 0;
    let mut query_span: u64 = 0;
    let mut ops = 0usize;
    for op in cigar.iter() {
        let op = match op {
            Ok(o) => o,
            Err(_) => break,
        };
        let len = op.len() as u64;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                ref_span = ref_span.saturating_add(len);
                query_span = query_span.saturating_add(len);
            }
            Kind::Insertion | Kind::SoftClip => {
                query_span = query_span.saturating_add(len);
            }
            Kind::Deletion | Kind::Skip => {
                ref_span = ref_span.saturating_add(len);
            }
            Kind::HardClip | Kind::Pad => {}
        }
        ops += 1;
        if ops > 4096 {
            break;
        }
    }

    // Sanity: both spans must remain inside u64 (saturating_add guarantees
    // this trivially, but the assertion documents the invariant).
    assert!(ref_span <= u64::MAX);
    assert!(query_span <= u64::MAX);
});
