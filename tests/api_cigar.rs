use rubam::api::Cigar;

#[test]
fn variant_lengths() {
    // Every variant must carry a u32 length, and matching is exhaustive.
    let cs = [
        Cigar::Match(5),
        Cigar::Ins(3),
        Cigar::Del(2),
        Cigar::RefSkip(100),
        Cigar::Equal(7),
        Cigar::Diff(1),
        Cigar::SoftClip(4),
        Cigar::HardClip(6),
        Cigar::Pad(0),
    ];
    for c in &cs {
        let len = match c {
            Cigar::Match(n) => *n,
            Cigar::Ins(n) => *n,
            Cigar::Del(n) => *n,
            Cigar::RefSkip(n) => *n,
            Cigar::Equal(n) => *n,
            Cigar::Diff(n) => *n,
            Cigar::SoftClip(n) => *n,
            Cigar::HardClip(n) => *n,
            Cigar::Pad(n) => *n,
        };
        // pure smoke-check: each value matches the constructor input
        let expected = match c {
            Cigar::Match(_) => 5,
            Cigar::Ins(_) => 3,
            Cigar::Del(_) => 2,
            Cigar::RefSkip(_) => 100,
            Cigar::Equal(_) => 7,
            Cigar::Diff(_) => 1,
            Cigar::SoftClip(_) => 4,
            Cigar::HardClip(_) => 6,
            Cigar::Pad(_) => 0,
        };
        assert_eq!(len, expected);
    }
}

#[test]
fn from_noodles_round_trip() {
    use noodles::sam::alignment::record::cigar::op::Kind;
    let pairs = [
        (Kind::Match, Cigar::Match(10)),
        (Kind::Insertion, Cigar::Ins(10)),
        (Kind::Deletion, Cigar::Del(10)),
        (Kind::Skip, Cigar::RefSkip(10)),
        (Kind::SequenceMatch, Cigar::Equal(10)),
        (Kind::SequenceMismatch, Cigar::Diff(10)),
        (Kind::SoftClip, Cigar::SoftClip(10)),
        (Kind::HardClip, Cigar::HardClip(10)),
        (Kind::Pad, Cigar::Pad(10)),
    ];
    for (kind, expected) in pairs {
        let c = Cigar::from_noodles_kind(kind, 10);
        assert_eq!(c, expected);
    }
}

#[test]
fn len_method_returns_u32() {
    assert_eq!(Cigar::Match(42).len(), 42);
    assert_eq!(Cigar::RefSkip(1_000_000).len(), 1_000_000);
}
