//! `rubam-depth` — pure-Rust CLI for per-base depth.
//!
//! Equivalent to `rubam depth` but without any Python in the hot path:
//! every byte of TSV output is produced by Rust, making this the fair
//! comparator against `samtools depth -a`.

use std::env;
use std::io::{self, BufWriter, Write};

use rubam::depth::compute_depths_native;

fn print_usage() {
    eprintln!(
        "rubam-depth — pure-Rust per-base depth\n\n\
         Usage: rubam-depth <bam> <chrom> <start> <end> \\\n\
         \t[-n THREADS] [-Q MIN_MAPQ] [-q MIN_BQ] [-d MAX_DEPTH] [-t STEP]\n\n\
         Output: chrom \\t pos \\t depth   (samtools-depth-compatible)"
    );
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.len() < 4 || args.iter().any(|a| a == "-h" || a == "--help") {
        print_usage();
        std::process::exit(if args.len() < 4 { 2 } else { 0 });
    }

    let bam = &args[0];
    let chrom = args[1].clone();
    let start: u64 = args[2].parse().unwrap_or_else(|_| die("invalid start"));
    let end: u64 = args[3].parse().unwrap_or_else(|_| die("invalid end"));

    let mut threads: usize = 12;
    let mut min_mapq: u8 = 0;
    let mut min_bq: u8 = 13;
    let mut max_depth: u32 = 8000;
    let mut step: u64 = 1;

    let mut i = 4;
    while i < args.len() {
        let next = || {
            args.get(i + 1)
                .cloned()
                .unwrap_or_else(|| die("missing value"))
        };
        match args[i].as_str() {
            "-n" | "--num-threads" => {
                threads = next().parse().unwrap_or(12);
                i += 2;
            }
            "-Q" | "--min-mapq" => {
                min_mapq = next().parse().unwrap_or(0);
                i += 2;
            }
            "-q" | "--min-bq" => {
                min_bq = next().parse().unwrap_or(13);
                i += 2;
            }
            "-d" | "--max-depth" => {
                max_depth = next().parse().unwrap_or(8000);
                i += 2;
            }
            "-t" | "--step" => {
                step = next().parse().unwrap_or(1);
                i += 2;
            }
            other => die(&format!("unknown option: {other}")),
        }
    }

    let (positions, depths) = compute_depths_native(
        bam, &chrom, start, end, step, min_mapq, min_bq, max_depth, threads,
    )
    .unwrap_or_else(|e| die(&e));

    let stdout = io::stdout();
    let mut out = BufWriter::with_capacity(1 << 20, stdout.lock());
    let mut buf = itoa_buf();
    for (p, d) in positions.iter().zip(depths.iter()) {
        out.write_all(chrom.as_bytes()).unwrap();
        out.write_all(b"\t").unwrap();
        out.write_all(write_u64(&mut buf, *p)).unwrap();
        out.write_all(b"\t").unwrap();
        out.write_all(write_u32(&mut buf, *d)).unwrap();
        out.write_all(b"\n").unwrap();
    }
    out.flush().unwrap();
}

fn die(msg: &str) -> ! {
    eprintln!("rubam-depth: {msg}");
    std::process::exit(2)
}

// Tiny integer formatter — avoids pulling in itoa crate.
fn itoa_buf() -> [u8; 32] {
    [0u8; 32]
}
fn write_u64(buf: &mut [u8; 32], mut n: u64) -> &[u8] {
    if n == 0 {
        buf[31] = b'0';
        return &buf[31..];
    }
    let mut i = 32;
    while n > 0 {
        i -= 1;
        buf[i] = b'0' + (n % 10) as u8;
        n /= 10;
    }
    &buf[i..]
}
fn write_u32(buf: &mut [u8; 32], n: u32) -> &[u8] {
    write_u64(buf, n as u64)
}
