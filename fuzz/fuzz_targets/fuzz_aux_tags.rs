#![no_main]
//! Fuzz target: arbitrary bytes → SAM record (with the aux/optional-fields
//! section coming from the fuzzer) → walk the aux iterator. This exercises
//! the same aux-tag decoder that `rubam::api::aux_data::aux_from_noodles`
//! consumes in `src/api/aux_data.rs`.
//!
//! We feed the fuzzer bytes as the optional-fields portion of a synthetic
//! SAM line. This is the cheapest way to drive the noodles aux parser
//! without depending on `pub(super)` constructors.

use libfuzzer_sys::fuzz_target;
use std::io::{BufReader, Cursor};

use noodles::sam;

/// Minimal SAM header sufficient to satisfy `read_header()`. One SQ line
/// referencing the contig we use in the synthetic record.
const HEADER: &[u8] = b"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n";

/// Mandatory fields of a SAM record up to (but not including) the
/// optional-fields portion. Trailing TAB separates from the aux block.
const RECORD_PREFIX: &[u8] = b"r1\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI\t";

fuzz_target!(|data: &[u8]| {
    // Bound the input length to keep iterations fast.
    if data.len() > 8192 {
        return;
    }
    // Reject embedded line terminators in the aux blob — those would
    // accidentally start a new record and obscure the surface we're
    // actually trying to exercise.
    if data.iter().any(|&b| b == b'\n' || b == b'\r') {
        return;
    }

    let mut buf = Vec::with_capacity(HEADER.len() + RECORD_PREFIX.len() + data.len() + 1);
    buf.extend_from_slice(HEADER);
    buf.extend_from_slice(RECORD_PREFIX);
    buf.extend_from_slice(data);
    buf.push(b'\n');

    let cursor = Cursor::new(buf);
    let mut reader = sam::io::Reader::new(BufReader::new(cursor));

    // If the header doesn't parse the fuzzer mutated something unexpected;
    // skip silently.
    let _header = match reader.read_header() {
        Ok(h) => h,
        Err(_) => return,
    };

    let mut record = sam::Record::default();
    match reader.read_record(&mut record) {
        Ok(0) | Err(_) => return,
        Ok(_) => {}
    }

    // Walk every (tag, value) pair the aux parser produces. Errors are
    // expected for malformed inputs; panics are the bugs we hunt.
    for entry in record.data().iter() {
        match entry {
            Ok((tag, value)) => {
                // Force the value to be touched: format the tag, and pattern
                // match through every Value variant the parser can emit.
                let _ = format!("{:?}", tag);
                use noodles::sam::alignment::record::data::field::Value as V;
                match value {
                    V::Character(c) => { let _ = c; }
                    V::Int8(n)   => { let _ = n; }
                    V::UInt8(n)  => { let _ = n; }
                    V::Int16(n)  => { let _ = n; }
                    V::UInt16(n) => { let _ = n; }
                    V::Int32(n)  => { let _ = n; }
                    V::UInt32(n) => { let _ = n; }
                    V::Float(n)  => { let _ = n; }
                    V::String(s) => { let _ = s.len(); }
                    V::Hex(s)    => { let _ = s.len(); }
                    V::Array(_a) => {
                        // Don't iterate Array contents here — that would
                        // require matching every Array subtype. The Value
                        // construction itself already crosses the parser.
                    }
                }
            }
            Err(_) => break,
        }
    }
});
