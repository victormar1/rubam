use rubam::api::AlignmentFile;

const FIXTURE: &str = "tests/example.bam";

#[test]
fn open_and_iter_count() {
    let mut bam = AlignmentFile::open(FIXTURE).unwrap();
    let n = bam.records().filter(|r| r.is_ok()).count();
    assert!(n > 0);
}

#[test]
fn header_has_chr1() {
    let bam = AlignmentFile::open(FIXTURE).unwrap();
    let h = bam.header();
    assert_eq!(h.tid2name(0), Some(&b"chr1"[..]));
}

#[test]
fn open_without_bai_works() {
    use std::fs;
    let tmp = std::env::temp_dir().join("rubam_no_bai.bam");
    fs::copy(FIXTURE, &tmp).unwrap();
    // Make sure no .bai is present
    let _ = fs::remove_file(format!("{}.bai", tmp.display()));
    let mut bam = AlignmentFile::open(&tmp).unwrap();
    let n = bam.records().filter(|r| r.is_ok()).count();
    assert!(n > 0);
    let _ = fs::remove_file(&tmp);
}

#[test]
fn open_missing_file_is_err() {
    let r = AlignmentFile::open("does_not_exist_zzz.bam");
    assert!(r.is_err());
}
