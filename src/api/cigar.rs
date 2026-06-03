//! Cigar — rust-htslib-compatible enum with the 9 SAM CIGAR ops.

use noodles::sam::alignment::record::cigar::op::Kind;

/// A single CIGAR operation. Variant names match `rust_htslib::bam::record::Cigar`
/// for drop-in compatibility with HARMOS-style code that pattern-matches:
///
/// ```ignore
/// match op {
///     Cigar::Match(n)    => ref_pos += n as u64,
///     Cigar::RefSkip(n)  => ref_pos += n as u64,
///     ...
/// }
/// ```
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Cigar {
    /// `M` — alignment match (consumes both ref and query).
    Match(u32),
    /// `I` — insertion in query (consumes query only).
    Ins(u32),
    /// `D` — deletion from reference (consumes reference only).
    Del(u32),
    /// `N` — skipped region from reference (consumes reference only).
    /// Used for spliced reads (RNA-seq introns).
    RefSkip(u32),
    /// `=` — sequence match (consumes both).
    Equal(u32),
    /// `X` — sequence mismatch (consumes both).
    Diff(u32),
    /// `S` — soft clip (consumes query only).
    SoftClip(u32),
    /// `H` — hard clip (consumes neither).
    HardClip(u32),
    /// `P` — padding (consumes neither).
    Pad(u32),
}

impl Cigar {
    /// The length (in bp) of this op.
    pub fn len(&self) -> u32 {
        match *self {
            Cigar::Match(n)
            | Cigar::Ins(n)
            | Cigar::Del(n)
            | Cigar::RefSkip(n)
            | Cigar::Equal(n)
            | Cigar::Diff(n)
            | Cigar::SoftClip(n)
            | Cigar::HardClip(n)
            | Cigar::Pad(n) => n,
        }
    }

    /// `true` iff this op is zero-length (rare but legal in BAM).
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Construct from a `noodles` `Kind` plus a length. Used by the bridge
    /// from `noodles_bam::Record::cigar()`.
    pub fn from_noodles_kind(kind: Kind, len: u32) -> Self {
        match kind {
            Kind::Match => Cigar::Match(len),
            Kind::Insertion => Cigar::Ins(len),
            Kind::Deletion => Cigar::Del(len),
            Kind::Skip => Cigar::RefSkip(len),
            Kind::SequenceMatch => Cigar::Equal(len),
            Kind::SequenceMismatch => Cigar::Diff(len),
            Kind::SoftClip => Cigar::SoftClip(len),
            Kind::HardClip => Cigar::HardClip(len),
            Kind::Pad => Cigar::Pad(len),
        }
    }

    /// `true` iff this op consumes the reference (M, D, N, =, X).
    pub fn consumes_reference(&self) -> bool {
        matches!(
            self,
            Cigar::Match(_) | Cigar::Del(_) | Cigar::RefSkip(_) | Cigar::Equal(_) | Cigar::Diff(_)
        )
    }

    /// `true` iff this op consumes the query (M, I, S, =, X).
    pub fn consumes_query(&self) -> bool {
        matches!(
            self,
            Cigar::Match(_) | Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::Equal(_) | Cigar::Diff(_)
        )
    }
}
