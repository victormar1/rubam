//! Cargo integration tests for the rubam-bcftools shadow CLI binary.
//!
//! Addresses Reviewer 2 M5 ('bcftools shadow CLI has zero cargo-test
//! coverage'). Each test invokes the binary as a subprocess and checks
//! exit code + a selective output substring.

use std::path::PathBuf;
use std::process::Command;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rubam-bcftools"))
}

fn fixture() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("tests/data/validation_3sample_100rec.vcf.gz");
    p
}

fn tmpdir(name: &str) -> PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!("rubam_bcftools_test_{}", name));
    let _ = std::fs::remove_dir_all(&p);
    std::fs::create_dir_all(&p).unwrap();
    p
}

#[test]
fn help_top_level() {
    let out = Command::new(binary()).output().unwrap();
    assert_eq!(out.status.code(), Some(2)); // no subcommand → 2
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("rubam-bcftools"));
    assert!(stderr.contains("view"));
    assert!(stderr.contains("norm"));
}

#[test]
fn help_explicit_flag() {
    let out = Command::new(binary()).arg("--help").output().unwrap();
    assert_eq!(out.status.code(), Some(0));
}

#[test]
fn unknown_subcommand_exits_2() {
    let out = Command::new(binary()).arg("frobnicate").output().unwrap();
    assert_eq!(out.status.code(), Some(2));
}

#[test]
fn view_count_matches_input() {
    let tmp = tmpdir("view_count");
    let out_path = tmp.join("out.vcf");
    let res = Command::new(binary())
        .args([
            "view",
            fixture().to_str().unwrap(),
            "-o",
            out_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert_eq!(
        res.status.code(),
        Some(0),
        "stderr: {}",
        String::from_utf8_lossy(&res.stderr)
    );
    assert!(out_path.exists(), "view did not write the output");
    let body = std::fs::read_to_string(&out_path).unwrap();
    let recs = body.lines().filter(|l| !l.starts_with('#')).count();
    assert_eq!(recs, 100, "expected 100 records, got {recs}");
}

#[test]
fn view_region_filter() {
    let tmp = tmpdir("view_region");
    let out_path = tmp.join("out.vcf");
    let res = Command::new(binary())
        .args([
            "view",
            "-r",
            "chr1",
            fixture().to_str().unwrap(),
            "-o",
            out_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert_eq!(res.status.code(), Some(0));
    let body = std::fs::read_to_string(&out_path).unwrap();
    let recs = body.lines().filter(|l| !l.starts_with('#')).count();
    // chr1, chr2, chr3 are present; chr1 should be a strict subset
    assert!(
        recs > 0 && recs < 100,
        "expected chr1-only subset, got {recs}"
    );
}

#[test]
fn query_format_string() {
    let tmp = tmpdir("query");
    let out_path = tmp.join("out.tsv");
    let res = Command::new(binary())
        .args([
            "query",
            "-f",
            "%CHROM\t%POS\n",
            "-o",
            out_path.to_str().unwrap(),
            fixture().to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert_eq!(res.status.code(), Some(0));
    let body = std::fs::read_to_string(&out_path).unwrap();
    let lines = body.lines().count();
    assert_eq!(lines, 100);
    // Each line should have one tab between CHROM and POS
    let first = body.lines().next().unwrap();
    assert!(first.contains('\t'));
}

#[test]
fn sort_round_trip_count() {
    let tmp = tmpdir("sort");
    let out_path = tmp.join("sorted.vcf");
    let res = Command::new(binary())
        .args([
            "sort",
            fixture().to_str().unwrap(),
            "-o",
            out_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert_eq!(res.status.code(), Some(0));
    let body = std::fs::read_to_string(&out_path).unwrap();
    let recs = body.lines().filter(|l| !l.starts_with('#')).count();
    assert_eq!(recs, 100);
}

#[test]
fn stats_emits_summary() {
    let res = Command::new(binary())
        .args(["stats", fixture().to_str().unwrap()])
        .output()
        .unwrap();
    assert_eq!(res.status.code(), Some(0));
    let stdout = String::from_utf8_lossy(&res.stdout);
    assert!(stdout.contains("total_records\t100"));
    assert!(stdout.contains("snps"));
    assert!(stdout.contains("ts_tv_ratio"));
}

#[test]
fn index_writes_tbi() {
    use std::fs;
    let tmp = tmpdir("index");
    let bgz = tmp.join("copy.vcf.gz");
    fs::copy(fixture(), &bgz).unwrap();
    let res = Command::new(binary())
        .args(["index", bgz.to_str().unwrap()])
        .output()
        .unwrap();
    assert_eq!(
        res.status.code(),
        Some(0),
        "stderr: {}",
        String::from_utf8_lossy(&res.stderr)
    );
    let tbi = bgz.with_extension("gz.tbi");
    assert!(tbi.exists(), "index did not produce a .tbi");
}

#[test]
fn missing_input_exits_nonzero() {
    let res = Command::new(binary())
        .args(["view", "/does/not/exist.vcf", "-o", "/tmp/dummy_rubam.vcf"])
        .output()
        .unwrap();
    assert_ne!(res.status.code(), Some(0));
}
