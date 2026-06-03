// tests/api_aligned_segment.rs
use rubam::api::{AlignedSegment, AlignmentFile, Cigar};

const FIXTURE: &str = "tests/example.bam";

fn first_record() -> AlignedSegment {
    let mut bam = AlignmentFile::open(FIXTURE).unwrap();
    bam.records().next().unwrap().unwrap()
}

#[test]
fn fields_basic() {
    let r = first_record();
    assert!(!r.qname().is_empty());
    assert!(r.tid() >= 0);
    assert!(r.pos() >= 0);
    // mapq() returns u8 so it is always <= 255 by type; the original assert
    // was a tautology and clippy rightly denies it.
    let _ = r.mapq();
    assert_eq!(r.seq().as_bytes().len(), r.seq_len());
    assert_eq!(r.qual().len(), r.seq_len());
}

#[test]
fn six_flag_accessors() {
    let r = first_record();
    let _ = r.is_unmapped();
    let _ = r.is_secondary();
    let _ = r.is_supplementary();
    let _ = r.is_duplicate();
    let _ = r.is_proper_pair();
    let _ = r.is_reverse();
    // example.bam has only proper-paired primary reads:
    assert!(r.is_proper_pair());
    assert!(!r.is_secondary());
    assert!(!r.is_supplementary());
}

#[test]
fn cigar_walks_in_enum_form() {
    let r = first_record();
    let ops: Vec<Cigar> = r.cigar().collect::<Result<_, _>>().unwrap();
    assert!(!ops.is_empty());
    // Total ref consumption + clips equals the read length on the query side
    // (only meaningful if the read is a clean M-only alignment, which is the
    // common case in example.bam).
    let query_consumed: u32 = ops
        .iter()
        .filter(|op| op.consumes_query())
        .map(|op| op.len())
        .sum();
    assert_eq!(query_consumed as usize, r.seq_len());
}

#[test]
fn qual_is_raw_phred_not_ascii_plus_33() {
    let r = first_record();
    let q = r.qual();
    if !q.is_empty() {
        // Raw phred values are in [0, 93]; ASCII +33 would put them in [33, 126].
        let max = *q.iter().max().unwrap();
        assert!(
            max <= 93,
            "qual seems to be ASCII +33 encoded (max={})",
            max
        );
    }
}
