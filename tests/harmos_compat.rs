// tests/harmos_compat.rs
//
// Integration tests asserting the rubam Rust API contract that HARMOS depends on.
// Drop-in compatible with the rust-htslib surface they migrate away from.
//
// NOTE: harmos_compat.sam is committed coordinate-sorted, so the BAM file order is:
//   [0] r1  flag=99    chr21:1000  (paired, proper, read1)
//   [1] r2  flag=147   chr21:1100  (paired, proper, read2, reverse)
//   [2] r4  flag=163   chr21:1900  (paired, proper, read2)  -- swapped vs SAM source
//   [3] r3  flag=83    chr21:2000  (paired, reverse, read2) -- SA tag here
//   [4] r5  flag=0     chr21:3000  (unpaired, complex cigar)
//   [5] r6  flag=1024  chr21:3500  (duplicate)
//   [6] r7  flag=256   chr21:3500  (secondary)
//   [7] r8  flag=2048  chr21:3500  (supplementary)
//   [8] r10 flag=16    chr22:5000  (reverse)
//   [9] r9  flag=77    *:0         (unmapped, mate unmapped)

use rubam::api::{AlignmentFile, Aux, Cigar};

const FIXTURE: &str = "tests/data/harmos_compat.bam";

#[test]
fn open_and_iter_count() {
    let mut bam = AlignmentFile::open(FIXTURE).expect("open");
    let n = bam.records().filter(|r| r.is_ok()).count();
    assert_eq!(n, 10, "fixture must yield exactly 10 records");
}

#[test]
fn header_tid_round_trip() {
    let bam = AlignmentFile::open(FIXTURE).expect("open");
    let h = bam.header();
    assert_eq!(h.target_count(), 2);
    assert_eq!(h.tid2name(0), Some(&b"chr21"[..]));
    assert_eq!(h.tid2name(1), Some(&b"chr22"[..]));
    assert_eq!(h.target_len(0), Some(46_709_983));
    assert_eq!(h.target_len(1), Some(50_818_468));
    assert_eq!(h.tid2name(99), None);
}

#[test]
fn record_fields_basic() {
    let mut bam = AlignmentFile::open(FIXTURE).expect("open");
    let r = bam.records().next().unwrap().unwrap();
    assert_eq!(r.qname(), b"r1");
    assert_eq!(r.tid(), 0);
    assert_eq!(r.pos(), 999); // 0-based
    assert_eq!(r.mapq(), 60);
    assert_eq!(r.seq_len(), 50);
    assert_eq!(r.seq().as_bytes().len(), 50);
    assert_eq!(r.qual().len(), 50);
}

#[test]
fn record_flags_all_six() {
    let mut bam = AlignmentFile::open(FIXTURE).expect("open");
    let records: Vec<_> = bam.records().filter_map(Result::ok).collect();
    // r1 at [0] (flag 99 = paired+proper+mate_reverse+read1)
    assert!(records[0].is_paired());
    assert!(records[0].is_proper_pair());
    assert!(!records[0].is_reverse());
    assert!(records[0].is_mate_reverse());
    // r6 at [5] (flag 1024 = duplicate)
    assert!(records[5].is_duplicate());
    // r7 at [6] (flag 256 = secondary)
    assert!(records[6].is_secondary());
    // r8 at [7] (flag 2048 = supplementary)
    assert!(records[7].is_supplementary());
    // r9 at [9] (flag 77 = paired+unmapped+mate_unmapped+read1)
    assert!(records[9].is_unmapped());
    assert!(records[9].is_mate_unmapped());
    // r10 at [8] (flag 16 = reverse only)
    assert!(records[8].is_reverse());
    assert!(!records[8].is_paired());
}

#[test]
fn cigar_all_nine_variants() {
    let mut bam = AlignmentFile::open(FIXTURE).expect("open");
    let records: Vec<_> = bam.records().filter_map(Result::ok).collect();
    // r5 at [4] has cigar 5S20M5I10M5N5=5X
    let r5 = &records[4];
    let ops: Vec<Cigar> = r5.cigar().collect::<Result<_, _>>().expect("cigar parse");
    assert_eq!(
        ops,
        vec![
            Cigar::SoftClip(5),
            Cigar::Match(20),
            Cigar::Ins(5),
            Cigar::Match(10),
            Cigar::RefSkip(5),
            Cigar::Equal(5),
            Cigar::Diff(5),
        ]
    );
}

#[test]
fn aux_sa_string() {
    let mut bam = AlignmentFile::open(FIXTURE).expect("open");
    let records: Vec<_> = bam.records().filter_map(Result::ok).collect();
    // r3 at [3] has SA:Z:chr22,5000,+,30M20S,60,0;
    let r3 = &records[3];
    match r3.aux(b"SA").expect("SA tag present") {
        Aux::String(s) => assert_eq!(s, "chr22,5000,+,30M20S,60,0;"),
        other => panic!("expected Aux::String, got {:?}", other),
    }
}

#[test]
fn aux_nm_int() {
    let mut bam = AlignmentFile::open(FIXTURE).expect("open");
    let r1 = bam.records().next().unwrap().unwrap();
    // r1 has NM:i:5 — could decode as I8/U8/I32/U32 depending on storage
    let nm = r1.aux(b"NM").expect("NM tag present");
    let v: i64 = match nm {
        Aux::I8(n) => n as i64,
        Aux::U8(n) => n as i64,
        Aux::I16(n) => n as i64,
        Aux::U16(n) => n as i64,
        Aux::I32(n) => n as i64,
        Aux::U32(n) => n as i64,
        other => panic!("expected integer Aux variant, got {:?}", other),
    };
    assert_eq!(v, 5);
}

#[test]
fn no_indexed_reader_path_works() {
    // HARMOS opens BAMs that may not have a .bai; rubam must succeed.
    use std::fs;
    let tmp = std::env::temp_dir().join("rubam_harmos_no_bai.bam");
    fs::copy(FIXTURE, &tmp).unwrap();
    let _ = fs::remove_file(format!("{}.bai", tmp.display()));
    let mut bam = AlignmentFile::open(&tmp).expect("open without .bai");
    let n = bam.records().filter(|r| r.is_ok()).count();
    assert_eq!(n, 10);
    let _ = fs::remove_file(&tmp);
}
