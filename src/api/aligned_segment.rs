//! AlignedSegment — single BAM record with HARMOS-shaped accessors.

use noodles::bam;

use super::cigar::Cigar;
use super::error::{Error, Result};

/// A single BAM alignment record. Constructed by `AlignmentFile::records()`.
///
/// Shape mirrors `rust_htslib::bam::Record` for HARMOS drop-in.
pub struct AlignedSegment {
    pub(crate) record: bam::Record,
    pub(crate) header: super::header::Header,
}

impl AlignedSegment {
    pub(crate) fn new(record: bam::Record, header: super::header::Header) -> Self {
        Self { record, header }
    }

    // ---------- identity / location ---------- //

    /// Read name (QNAME). Empty bytes when no name is present.
    pub fn qname(&self) -> &[u8] {
        self.record.name().map(|n| n.as_ref()).unwrap_or(b"")
    }

    /// Reference sequence ID (TID). Returns `-1` for unmapped reads or
    /// records without a reference sequence.
    pub fn tid(&self) -> i32 {
        match self.record.reference_sequence_id() {
            Some(Ok(id)) => id as i32,
            _ => -1,
        }
    }

    /// 0-based start position. Returns `-1` for unmapped reads.
    pub fn pos(&self) -> i64 {
        match self.record.alignment_start() {
            Some(Ok(p)) => (p.get() as i64) - 1,
            _ => -1,
        }
    }

    /// Mapping quality (0-255). 255 means missing.
    pub fn mapq(&self) -> u8 {
        self.record
            .mapping_quality()
            .map(|q| q.get())
            .unwrap_or(255)
    }

    /// Length of the read sequence.
    pub fn seq_len(&self) -> usize {
        self.record.sequence().len()
    }

    /// Sequence as a byte buffer (`&[u8]` of `b'A'`/`b'C'`/`b'G'`/`b'T'`/`b'N'`).
    pub fn seq(&self) -> SeqBytes<'_> {
        SeqBytes {
            record: &self.record,
        }
    }

    /// Per-base quality scores in **raw phred** (0-93), NOT ASCII +33.
    pub fn qual(&self) -> Vec<u8> {
        self.record.quality_scores().iter().collect()
    }

    // ---------- flags ---------- //

    fn flag_bit(&self, mask: u16) -> bool {
        self.record.flags().bits() & mask != 0
    }

    pub fn is_paired(&self) -> bool {
        self.flag_bit(0x1)
    }
    pub fn is_proper_pair(&self) -> bool {
        self.flag_bit(0x2)
    }
    pub fn is_unmapped(&self) -> bool {
        self.flag_bit(0x4)
    }
    pub fn is_mate_unmapped(&self) -> bool {
        self.flag_bit(0x8)
    }
    pub fn is_reverse(&self) -> bool {
        self.flag_bit(0x10)
    }
    pub fn is_mate_reverse(&self) -> bool {
        self.flag_bit(0x20)
    }
    pub fn is_read1(&self) -> bool {
        self.flag_bit(0x40)
    }
    pub fn is_read2(&self) -> bool {
        self.flag_bit(0x80)
    }
    pub fn is_secondary(&self) -> bool {
        self.flag_bit(0x100)
    }
    pub fn is_qcfail(&self) -> bool {
        self.flag_bit(0x200)
    }
    pub fn is_duplicate(&self) -> bool {
        self.flag_bit(0x400)
    }
    pub fn is_supplementary(&self) -> bool {
        self.flag_bit(0x800)
    }

    // ---------- CIGAR ---------- //

    /// Iterator over CIGAR ops as `Cigar` enum variants (HARMOS-style match-
    /// friendly). Yields `Result` because CIGAR decoding can fail on
    /// malformed BAMs.
    ///
    /// Note: ops are decoded eagerly into a `Vec` to avoid lifetime issues
    /// with noodles' Cigar temporary under Rust 2024 capture rules.
    pub fn cigar(&self) -> impl Iterator<Item = Result<Cigar>> + '_ {
        let ops: Vec<Result<Cigar>> = self
            .record
            .cigar()
            .iter()
            .map(|res| {
                res.map_err(|e| Error::Cigar(e.to_string()))
                    .map(|op| Cigar::from_noodles_kind(op.kind(), op.len() as u32))
            })
            .collect();
        ops.into_iter()
    }

    // ---------- aux tags ---------- //

    /// Look up an auxiliary tag by 2-byte name. Returns `Aux::String` /
    /// `Aux::I32` / etc. depending on the encoded type.
    pub fn aux<'r>(
        &'r self,
        tag: &[u8],
    ) -> std::result::Result<super::aux_data::Aux<'r>, super::aux_data::AuxError> {
        use super::aux_data::AuxError;
        if tag.len() != 2 {
            return Err(AuxError::BadTagLength(tag.len()));
        }
        let target: [u8; 2] = [tag[0], tag[1]];
        let data = self.record.data();
        for entry in data.iter() {
            let (t, v) = entry.map_err(|e| {
                AuxError::Parse(String::from_utf8_lossy(&target).into_owned(), e.to_string())
            })?;
            if t.as_ref() == &target {
                // v0.3.2 Wave 3: no more Box::leak. Lifetimes flow purely from
                // the noodles value; the arena parameter is gone. Array tags
                // still return Unsupported until v0.4 introduces an owned arena.
                return super::aux_data::aux_from_noodles(target, v);
            }
        }
        Err(AuxError::NotFound(
            String::from_utf8_lossy(&target).into_owned(),
        ))
    }

    /// The header this record was read against.
    pub fn header(&self) -> &super::header::Header {
        &self.header
    }
}

/// Borrowed view over the read's sequence bytes. `as_bytes()` returns a
/// freshly allocated `Vec<u8>`; `iter()` yields each base as ASCII without
/// allocating.
pub struct SeqBytes<'a> {
    record: &'a bam::Record,
}

impl<'a> SeqBytes<'a> {
    /// Materialize the sequence as a `Vec<u8>` of ASCII bases (`A`/`C`/`G`/`T`/`N`).
    pub fn as_bytes(&self) -> Vec<u8> {
        let seq = self.record.sequence();
        let n = seq.len();
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            out.push(seq.get(i).unwrap_or(b'N'));
        }
        out
    }

    pub fn len(&self) -> usize {
        self.record.sequence().len()
    }

    pub fn is_empty(&self) -> bool {
        self.record.sequence().len() == 0
    }
}
