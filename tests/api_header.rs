use rubam::api::Header;

const FIXTURE: &str = "tests/example.bam";

fn open_header() -> Header {
    rubam::api::AlignmentFile::open(FIXTURE)
        .unwrap()
        .header()
        .clone()
}

#[test]
fn target_count_matches_references() {
    let h = open_header();
    assert_eq!(h.target_count() as usize, h.target_names().count());
}

#[test]
fn tid2name_chr1() {
    let h = open_header();
    assert_eq!(h.tid2name(0), Some(&b"chr1"[..]));
}

#[test]
fn target_len_chr1_positive() {
    let h = open_header();
    assert!(h.target_len(0).map(|l| l > 0).unwrap_or(false));
}

#[test]
fn tid2name_oob_is_none() {
    let h = open_header();
    assert_eq!(h.tid2name(99_999), None);
    assert_eq!(h.target_len(99_999), None);
}

#[test]
fn target_names_iterates_in_order() {
    let h = open_header();
    let names: Vec<&[u8]> = h.target_names().collect();
    assert!(!names.is_empty());
    assert_eq!(names[0], b"chr1");
}
