#![no_main]
//! Fuzz target: feed arbitrary bytes to `noodles::bam::io::Reader` and
//! ensure the parser does not panic. Errors are expected — only panics
//! and aborts count as findings.
//!
//! Coverage: BGZF wrapper + BAM header parser + BAM record decoder, i.e.
//! the exact entry points used by `rubam::api::AlignmentFile::open`.

use libfuzzer_sys::fuzz_target;
use std::io::Cursor;

use noodles::bam;

fuzz_target!(|data: &[u8]| {
    // bam::io::Reader::new wraps the inner reader in a bgzf::io::Reader
    // internally, so feeding raw arbitrary bytes exercises both the BGZF
    // block decoder and the BAM record parser.
    let cursor = Cursor::new(data);
    let mut reader = bam::io::Reader::new(cursor);

    // Header parse — most malformed inputs die here, which is fine: the
    // contract is "return Err, don't panic".
    let _ = match reader.read_header() {
        Ok(h) => h,
        Err(_) => return,
    };

    // If the header parsed, try to read up to a small bounded number of
    // records so we exercise the record decoder without letting the fuzzer
    // sit forever on a giant valid prefix.
    let mut record = bam::Record::default();
    for _ in 0..64 {
        match reader.read_record(&mut record) {
            Ok(0) => break,
            Ok(_) => {
                // Touch a few accessors so the fuzzer sees coverage on the
                // lazy-decode paths (CIGAR, sequence, quality_scores, data).
                let _ = record.flags();
                let _ = record.reference_sequence_id();
                let _ = record.alignment_start();
                let _ = record.cigar().len();
                let _ = record.sequence().len();
                let _ = record.quality_scores().len();
                let _ = record.data().iter().count();
            }
            Err(_) => break,
        }
    }
});
