//! `rubam samtools <subcommand>` shadow CLI.
//!
//! Mirrors a subset of the samtools CLI for drop-in replacement. Each
//! subcommand dispatches to a function in `rubam::tools::<subcmd>`.
//! The full subcommand wiring lands in tasks C2..C5; this skeleton
//! just sets up the dispatcher and `--help`.

use std::env;

fn print_top_help() {
    eprintln!(
        "rubam-samtools — drop-in samtools subcommand dispatcher\n\n\
         Usage: rubam-samtools <subcommand> [args...]\n\n\
         Subcommands:\n  \
           sort      Coordinate-sort a BAM\n  \
           index     Write BAI for a sorted BAM\n  \
           view      Region/flag/MAPQ filter and BAM output\n  \
           merge     Merge multiple sorted BAMs into one\n  \
           flagstat  samtools flagstat\n  \
           idxstats  samtools idxstats\n  \
           calmd     Recompute NM (and MD in v0.2.x)\n  \
           faidx     Build FASTA index / extract subsequence\n"
    );
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.is_empty() || args[0] == "-h" || args[0] == "--help" {
        print_top_help();
        std::process::exit(if args.is_empty() { 2 } else { 0 });
    }
    let sub = &args[0];
    let rest = &args[1..];
    let rc = match sub.as_str() {
        "sort" => sub_sort(rest),
        "index" => sub_index(rest),
        "view" => sub_view(rest),
        "merge" => sub_merge(rest),
        "flagstat" => sub_flagstat(rest),
        "idxstats" => sub_idxstats(rest),
        "calmd" => sub_calmd(rest),
        "faidx" => sub_faidx(rest),
        other => {
            eprintln!("rubam-samtools: unknown subcommand {other:?}");
            print_top_help();
            2
        }
    };
    std::process::exit(rc);
}

// Each sub_* below parses its own argv and calls the corresponding
// rubam::tools::* function. Filled in by tasks C2-C5.

fn sub_sort(args: &[String]) -> i32 {
    // Accepted forms:  sort [-o OUT] [-@ N] IN
    let mut output: Option<String> = None;
    let mut threads: usize = 1;
    let mut input: Option<String> = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-o" => {
                output = args.get(i + 1).cloned();
                if output.is_none() {
                    eprintln!("sort: -o needs an argument");
                    return 2;
                }
                i += 2;
            }
            "-@" => {
                threads = args.get(i + 1).and_then(|s| s.parse().ok()).unwrap_or(1);
                i += 2;
            }
            x if !x.starts_with('-') => {
                input = Some(x.to_string());
                i += 1;
            }
            other => {
                eprintln!("sort: unknown arg {other}");
                return 2;
            }
        }
    }
    let _ = threads; // accepted for samtools parity, ignored in v0.2
    let Some(input) = input else {
        eprintln!("sort: missing INPUT");
        return 2;
    };
    let Some(output) = output else {
        eprintln!("sort: missing -o OUTPUT");
        return 2;
    };
    if let Err(e) = rubam::tools::sort::sort_native(&input, &output) {
        eprintln!("sort: {e}");
        return 1;
    }
    0
}

fn sub_index(args: &[String]) -> i32 {
    let mut input: Option<String> = None;
    let mut csi = false;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-c" | "--csi" => {
                csi = true;
                i += 1;
            }
            x if !x.starts_with('-') => {
                input = Some(x.to_string());
                i += 1;
            }
            other => {
                eprintln!("index: unknown arg {other}");
                return 2;
            }
        }
    }
    if csi {
        eprintln!("index: CSI lands in v0.2.x; use BAI by omitting -c");
        return 2;
    }
    let Some(input) = input else {
        eprintln!("index: missing INPUT");
        return 2;
    };
    if let Err(e) = rubam::tools::index::index_native(&input) {
        eprintln!("index: {e}");
        return 1;
    }
    0
}
fn sub_view(args: &[String]) -> i32 {
    let mut input: Option<String> = None;
    let mut output: Option<String> = None;
    let mut region: Option<String> = None;
    let mut count_only = false;
    let mut want_bam_out = false;
    let mut min_mapq: u8 = 0;
    let mut flag_required: u16 = 0;
    let mut flag_filtered: u16 = 0;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-c" => {
                count_only = true;
                i += 1;
            }
            "-b" => {
                want_bam_out = true;
                i += 1;
            }
            "-o" => {
                output = args.get(i + 1).cloned();
                if output.is_none() {
                    eprintln!("view: -o needs an argument");
                    return 2;
                }
                i += 2;
            }
            "-q" => {
                min_mapq = args.get(i + 1).and_then(|s| s.parse().ok()).unwrap_or(0);
                i += 2;
            }
            "-f" => {
                flag_required = args.get(i + 1).and_then(|s| parse_flag(s)).unwrap_or(0);
                i += 2;
            }
            "-F" => {
                flag_filtered = args.get(i + 1).and_then(|s| parse_flag(s)).unwrap_or(0);
                i += 2;
            }
            x if !x.starts_with('-') => {
                if input.is_none() {
                    input = Some(x.to_string());
                } else {
                    region = Some(x.to_string());
                }
                i += 1;
            }
            other => {
                eprintln!("view: unknown arg {other}");
                return 2;
            }
        }
    }
    let _ = want_bam_out; // we always emit BAM when output is given.
    let Some(input) = input else {
        eprintln!("view: missing INPUT");
        return 2;
    };
    match rubam::tools::view::view_native(
        &input,
        region.as_deref(),
        output.as_deref(),
        min_mapq,
        flag_required,
        flag_filtered,
        count_only,
    ) {
        Ok(n) => {
            if count_only {
                println!("{n}");
            }
            0
        }
        Err(e) => {
            eprintln!("view: {e}");
            1
        }
    }
}

/// Parse a SAM flag value: accepts decimal or 0x-prefixed hex.
fn parse_flag(s: &str) -> Option<u16> {
    if let Some(stripped) = s.strip_prefix("0x").or_else(|| s.strip_prefix("0X")) {
        u16::from_str_radix(stripped, 16).ok()
    } else {
        s.parse().ok()
    }
}

fn sub_merge(args: &[String]) -> i32 {
    if args.len() < 2 {
        eprintln!("merge: usage: merge OUTPUT INPUT [INPUT ...]");
        return 2;
    }
    let output = args[0].clone();
    let inputs: Vec<&str> = args[1..].iter().map(|s| s.as_str()).collect();
    if let Err(e) = rubam::tools::merge::merge_native(&inputs, &output, true) {
        eprintln!("merge: {e}");
        return 1;
    }
    0
}
fn sub_flagstat(args: &[String]) -> i32 {
    let Some(input) = args.first() else {
        eprintln!("flagstat: missing INPUT");
        return 2;
    };
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    if let Err(e) = rubam::tools::flagstat::flagstat_native(input, &mut out) {
        eprintln!("flagstat: {e}");
        return 1;
    }
    0
}

fn sub_idxstats(args: &[String]) -> i32 {
    let Some(input) = args.first() else {
        eprintln!("idxstats: missing INPUT");
        return 2;
    };
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    if let Err(e) = rubam::tools::idxstats::idxstats_native(input, &mut out) {
        eprintln!("idxstats: {e}");
        return 1;
    }
    0
}
fn sub_calmd(args: &[String]) -> i32 {
    // Accepted: calmd [-b] INPUT REFERENCE — output goes to stdout (BAM).
    let mut input: Option<String> = None;
    let mut reference: Option<String> = None;
    let mut emit_bam = false;
    for a in args {
        match a.as_str() {
            "-b" => emit_bam = true,
            x if !x.starts_with('-') => {
                if input.is_none() {
                    input = Some(x.to_string());
                } else if reference.is_none() {
                    reference = Some(x.to_string());
                } else {
                    eprintln!("calmd: unexpected positional arg {x:?}");
                    return 2;
                }
            }
            other => {
                eprintln!("calmd: unknown arg {other}");
                return 2;
            }
        }
    }
    let _ = emit_bam; // we always emit BAM in v0.2
    let Some(input) = input else {
        eprintln!("calmd: missing INPUT");
        return 2;
    };
    let Some(reference) = reference else {
        eprintln!("calmd: missing REFERENCE");
        return 2;
    };

    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    let mut writer = noodles::bam::io::Writer::new(noodles::bgzf::io::Writer::new(&mut out));
    if let Err(e) = rubam::tools::calmd::calmd_native_to(&input, &reference, &mut writer) {
        eprintln!("calmd: {e}");
        return 1;
    }
    0
}

fn sub_faidx(args: &[String]) -> i32 {
    // Accepted: faidx FASTA [REGION ...]
    let Some(fasta) = args.first() else {
        eprintln!("faidx: missing FASTA");
        return 2;
    };
    let regions = &args[1..];
    if regions.is_empty() {
        if let Err(e) = rubam::tools::faidx::faidx_index_only(fasta) {
            eprintln!("faidx: {e}");
            return 1;
        }
        return 0;
    }
    for r in regions {
        match rubam::tools::faidx::faidx_subseq(fasta, r) {
            Ok((header, seq)) => {
                println!(">{header}");
                println!("{seq}");
            }
            Err(e) => {
                eprintln!("faidx: {e}");
                return 1;
            }
        }
    }
    0
}
