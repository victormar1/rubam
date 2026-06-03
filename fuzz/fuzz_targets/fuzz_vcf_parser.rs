#![no_main]
//! Fuzz target: feed arbitrary bytes to `noodles::vcf::io::Reader` and
//! ensure the parser does not panic. The reader is constructed from a
//! `BufReader` over the arbitrary bytes, matching how
//! `rubam::variant::VariantFile` opens plain VCF inputs.

use libfuzzer_sys::fuzz_target;
use std::io::{BufReader, Cursor};

use noodles::vcf;
// Bring the field-iterator traits into scope so we can call .iter()
// on AlternateBases / Ids / Filters / Samples returned by Record.
use noodles::vcf::variant::record::{
    AlternateBases as _, Filters as _, Ids as _,
};

fuzz_target!(|data: &[u8]| {
    let buf = BufReader::new(Cursor::new(data));
    let mut reader = vcf::io::Reader::new(buf);

    // VCF header: text-based, riddled with edge cases (key=value, INFO
    // definitions, contig lines, escaping). This is the highest-value
    // surface to fuzz.
    let header = match reader.read_header() {
        Ok(h) => h,
        Err(_) => return,
    };

    // If the header parsed, walk a bounded number of records using the
    // typed `Record` reader. Each call lex-parses one line and exposes
    // the field accessors below, which lazily parse on demand.
    let mut record = vcf::Record::default();
    for _ in 0..128 {
        match reader.read_record(&mut record) {
            Ok(0) => break,
            Ok(_) => {
                // Touch the lazy-parse accessors so the fuzzer gets
                // coverage on the per-field parsers (POS, REF, ALT, FILTER,
                // INFO, FORMAT/samples).
                let _ = record.reference_sequence_name();
                let _ = record.variant_start();
                let _ = record.reference_bases();
                let _ = record.alternate_bases().iter().count();
                let _ = record.ids().iter().count();
                let _ = record.filters().iter(&header).count();
                let _ = record.info().iter(&header).count();
                let _ = record.samples().iter().count();
            }
            Err(_) => break,
        }
    }
});
