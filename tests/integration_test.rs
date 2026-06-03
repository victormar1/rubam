//! Rust-side integration tests. The pyo3-driven path is tested from Python
//! (see `tests/test_core.py`); here we only make sure the BAM fixture is
//! readable through the noodles backend on the current target.

use noodles::bam;

#[test]
fn reads_example_bam_with_noodles() {
    let bam_path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/example.bam");
    let mut reader = bam::io::reader::Builder
        .build_from_path(bam_path)
        .expect("open BAM");
    let _header = reader.read_header().expect("read header");

    let mut count = 0usize;
    for result in reader.records() {
        let record = result.expect("read record");
        // Touch one method to make sure the record is decoded.
        let _ = record.flags();
        count += 1;
    }
    assert!(count > 0, "no records found in BAM file");
}
