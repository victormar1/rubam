//! `rubam bcftools <subcommand>` shadow CLI.
//!
//! Mirrors a subset of the bcftools CLI for drop-in replacement. Each
//! subcommand dispatches to a function in `rubam::tools::bcftools::<subcmd>`.
//! The 7 subcommands wired here cover Phase B of `paper/PLAN_v0.3.md`.

use std::env;
use std::path::Path;

use rubam::tools::bcftools;

fn print_top_help() {
    eprintln!(
        "rubam-bcftools — drop-in bcftools subcommand dispatcher\n\n\
         Usage: rubam-bcftools <subcommand> [args...]\n\n\
         Subcommands:\n  \
           view      VCF/BCF region/sample/format filter\n  \
           norm      Split multi-allelic + left-align indels\n  \
           concat    Glue sorted VCFs/BCFs together\n  \
           query     Format-string field extraction\n  \
           index     Build TBI / CSI index\n  \
           sort      Coordinate-sort VCF/BCF\n  \
           stats     Per-sample summary (Ts/Tv, het/hom, totals)\n"
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
        "view" => sub_view(rest),
        "norm" => sub_norm(rest),
        "concat" => sub_concat(rest),
        "query" => sub_query(rest),
        "index" => sub_index(rest),
        "sort" => sub_sort(rest),
        "stats" => sub_stats(rest),
        other => {
            eprintln!("rubam-bcftools: unknown subcommand {other:?}");
            print_top_help();
            2
        }
    };
    std::process::exit(rc);
}

/// Tag-based representation of -O so the CLI is uncoupled from any tool's
/// concrete `OutputFormat` enum (each subcommand defines its own copy).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FormatTag {
    Vcf,
    VcfGz,
    Bcf,
}

fn parse_format(s: &str) -> Result<FormatTag, i32> {
    match s {
        "v" | "vcf" => Ok(FormatTag::Vcf),
        "z" | "vcf.gz" => Ok(FormatTag::VcfGz),
        "b" | "bcf" => Ok(FormatTag::Bcf),
        _ => {
            eprintln!("unknown -O output type {s:?} (use v|z|b)");
            Err(2)
        }
    }
}

fn to_view_fmt(t: FormatTag) -> bcftools::view::OutputFormat {
    match t {
        FormatTag::Vcf => bcftools::view::OutputFormat::Vcf,
        FormatTag::VcfGz => bcftools::view::OutputFormat::VcfGz,
        FormatTag::Bcf => bcftools::view::OutputFormat::Bcf,
    }
}

fn to_sort_fmt(t: FormatTag) -> bcftools::sort::OutputFormat {
    match t {
        FormatTag::Vcf => bcftools::sort::OutputFormat::Vcf,
        FormatTag::VcfGz => bcftools::sort::OutputFormat::VcfGz,
        FormatTag::Bcf => bcftools::sort::OutputFormat::Bcf,
    }
}

fn to_concat_fmt(t: FormatTag) -> bcftools::concat::OutputFormat {
    match t {
        FormatTag::Vcf => bcftools::concat::OutputFormat::Vcf,
        FormatTag::VcfGz => bcftools::concat::OutputFormat::VcfGz,
        FormatTag::Bcf => bcftools::concat::OutputFormat::Bcf,
    }
}

fn to_norm_fmt(t: FormatTag) -> bcftools::norm::OutputFormat {
    match t {
        FormatTag::Vcf => bcftools::norm::OutputFormat::Vcf,
        FormatTag::VcfGz => bcftools::norm::OutputFormat::VcfGz,
        FormatTag::Bcf => bcftools::norm::OutputFormat::Bcf,
    }
}

fn sub_view(args: &[String]) -> i32 {
    let mut input: Option<String> = None;
    let mut region: Option<String> = None;
    let mut samples: Option<Vec<String>> = None;
    let mut output: Option<String> = None;
    let mut format_str = String::from("v");
    let mut header_only = false;
    let mut no_header = false;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-r" | "--region" => {
                region = args.get(i + 1).cloned();
                if region.is_none() {
                    eprintln!("view: -r needs an argument");
                    return 2;
                }
                i += 2;
            }
            "-s" | "--samples" => {
                let s = args.get(i + 1).cloned();
                if s.is_none() {
                    eprintln!("view: -s needs an argument");
                    return 2;
                }
                samples = Some(s.unwrap().split(',').map(String::from).collect());
                i += 2;
            }
            "-o" | "--output" => {
                output = args.get(i + 1).cloned();
                if output.is_none() {
                    eprintln!("view: -o needs an argument");
                    return 2;
                }
                i += 2;
            }
            "-O" | "--output-type" => {
                format_str = args.get(i + 1).cloned().unwrap_or_default();
                i += 2;
            }
            "-h" | "--header-only" => {
                header_only = true;
                i += 1;
            }
            "-H" | "--no-header" => {
                no_header = true;
                i += 1;
            }
            x if !x.starts_with('-') => {
                input = Some(x.to_string());
                i += 1;
            }
            other => {
                eprintln!("view: unknown arg {other}");
                return 2;
            }
        }
    }
    let Some(input) = input else {
        eprintln!("view: missing INPUT");
        return 2;
    };
    let Some(output) = output else {
        eprintln!("view: -o OUTPUT is required in v0.3");
        return 2;
    };
    let format = match parse_format(&format_str) {
        Ok(f) => to_view_fmt(f),
        Err(rc) => return rc,
    };
    let samples_slice = samples.as_deref();
    let out_path = Path::new(&output);
    match bcftools::view::view_native(
        &input,
        region.as_deref(),
        samples_slice,
        Some(out_path),
        format,
        header_only,
        no_header,
    ) {
        Ok(_) => 0,
        Err(e) => {
            eprintln!("view: {e}");
            1
        }
    }
}

fn sub_norm(args: &[String]) -> i32 {
    let mut input: Option<String> = None;
    let mut output: Option<String> = None;
    let mut format_str = String::from("v");
    let mut multiallelic = String::new();
    let mut reference: Option<String> = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-o" | "--output" => {
                output = args.get(i + 1).cloned();
                if output.is_none() {
                    eprintln!("norm: -o needs an argument");
                    return 2;
                }
                i += 2;
            }
            "-O" | "--output-type" => {
                format_str = args.get(i + 1).cloned().unwrap_or_default();
                i += 2;
            }
            "-m" => {
                multiallelic = args.get(i + 1).cloned().unwrap_or_default();
                i += 2;
            }
            "-f" | "--reference" => {
                reference = args.get(i + 1).cloned();
                if reference.is_none() {
                    eprintln!("norm: -f needs a path");
                    return 2;
                }
                i += 2;
            }
            x if !x.starts_with('-') => {
                input = Some(x.to_string());
                i += 1;
            }
            other => {
                eprintln!("norm: unknown arg {other}");
                return 2;
            }
        }
    }
    let Some(input) = input else {
        eprintln!("norm: missing INPUT");
        return 2;
    };
    let Some(output) = output else {
        eprintln!("norm: -o OUTPUT is required");
        return 2;
    };
    let format = match parse_format(&format_str) {
        Ok(f) => to_norm_fmt(f),
        Err(rc) => return rc,
    };
    let split = match multiallelic.as_str() {
        "" | "-" => multiallelic == "-",
        "+" => {
            eprintln!("norm: -m + (join) lands in v0.3.x");
            return 2;
        }
        other => {
            eprintln!("norm: unknown -m {other:?}");
            return 2;
        }
    };
    let ref_path = reference.as_deref().map(Path::new);
    match bcftools::norm::norm_native(&input, Path::new(&output), format, split, ref_path) {
        Ok((in_n, out_n, la_n)) => {
            eprintln!("norm: {in_n} in -> {out_n} out, {la_n} left-aligned");
            0
        }
        Err(e) => {
            eprintln!("norm: {e}");
            1
        }
    }
}

fn sub_concat(args: &[String]) -> i32 {
    let mut inputs: Vec<String> = Vec::new();
    let mut output: Option<String> = None;
    let mut format_str = String::from("v");
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-o" | "--output" => {
                output = args.get(i + 1).cloned();
                if output.is_none() {
                    eprintln!("concat: -o needs an argument");
                    return 2;
                }
                i += 2;
            }
            "-O" | "--output-type" => {
                format_str = args.get(i + 1).cloned().unwrap_or_default();
                i += 2;
            }
            x if !x.starts_with('-') => {
                inputs.push(x.to_string());
                i += 1;
            }
            other => {
                eprintln!("concat: unknown arg {other}");
                return 2;
            }
        }
    }
    if inputs.len() < 2 {
        eprintln!("concat: needs at least 2 inputs");
        return 2;
    }
    let Some(output) = output else {
        eprintln!("concat: -o OUTPUT is required");
        return 2;
    };
    let format = match parse_format(&format_str) {
        Ok(f) => to_concat_fmt(f),
        Err(rc) => return rc,
    };
    match bcftools::concat::concat_native(&inputs, Path::new(&output), format) {
        Ok(n) => {
            eprintln!("concat: wrote {n} records");
            0
        }
        Err(e) => {
            eprintln!("concat: {e}");
            1
        }
    }
}

fn sub_query(args: &[String]) -> i32 {
    let mut input: Option<String> = None;
    let mut format: Option<String> = None;
    let mut output: Option<String> = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-f" | "--format" => {
                format = args.get(i + 1).cloned();
                if format.is_none() {
                    eprintln!("query: -f needs an argument");
                    return 2;
                }
                i += 2;
            }
            "-o" | "--output" => {
                output = args.get(i + 1).cloned();
                if output.is_none() {
                    eprintln!("query: -o needs an argument");
                    return 2;
                }
                i += 2;
            }
            x if !x.starts_with('-') => {
                input = Some(x.to_string());
                i += 1;
            }
            other => {
                eprintln!("query: unknown arg {other}");
                return 2;
            }
        }
    }
    let Some(input) = input else {
        eprintln!("query: missing INPUT");
        return 2;
    };
    let Some(format) = format else {
        eprintln!("query: -f FORMAT is required");
        return 2;
    };
    let Some(output) = output else {
        eprintln!("query: -o OUTPUT is required in v0.3");
        return 2;
    };
    match bcftools::query::query_native(&input, &format, Path::new(&output)) {
        Ok(n) => {
            eprintln!("query: emitted {n} records");
            0
        }
        Err(e) => {
            eprintln!("query: {e}");
            1
        }
    }
}

fn sub_index(args: &[String]) -> i32 {
    let mut input: Option<String> = None;
    let mut csi = false;
    let mut force = false;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-c" | "--csi" => {
                csi = true;
                i += 1;
            }
            "-t" | "--tbi" => {
                csi = false;
                i += 1;
            }
            "-f" | "--force" => {
                force = true;
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
    let Some(input) = input else {
        eprintln!("index: missing INPUT");
        return 2;
    };
    let kind = if csi {
        bcftools::index::IndexKind::Csi
    } else {
        bcftools::index::IndexKind::Tbi
    };
    match bcftools::index::index_native(&input, kind, force) {
        Ok(p) => {
            eprintln!("index: wrote {}", p.display());
            0
        }
        Err(e) => {
            eprintln!("index: {e}");
            1
        }
    }
}

fn sub_sort(args: &[String]) -> i32 {
    let mut input: Option<String> = None;
    let mut output: Option<String> = None;
    let mut format_str = String::from("v");
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-o" | "--output" => {
                output = args.get(i + 1).cloned();
                if output.is_none() {
                    eprintln!("sort: -o needs an argument");
                    return 2;
                }
                i += 2;
            }
            "-O" | "--output-type" => {
                format_str = args.get(i + 1).cloned().unwrap_or_default();
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
    let Some(input) = input else {
        eprintln!("sort: missing INPUT");
        return 2;
    };
    let Some(output) = output else {
        eprintln!("sort: -o OUTPUT is required");
        return 2;
    };
    let format = match parse_format(&format_str) {
        Ok(f) => to_sort_fmt(f),
        Err(rc) => return rc,
    };
    match bcftools::sort::sort_native(&input, Path::new(&output), format) {
        Ok(n) => {
            eprintln!("sort: wrote {n} records");
            0
        }
        Err(e) => {
            eprintln!("sort: {e}");
            1
        }
    }
}

fn sub_stats(args: &[String]) -> i32 {
    let mut input: Option<String> = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            x if !x.starts_with('-') => {
                input = Some(x.to_string());
                i += 1;
            }
            other => {
                eprintln!("stats: unknown arg {other}");
                return 2;
            }
        }
    }
    let Some(input) = input else {
        eprintln!("stats: missing INPUT");
        return 2;
    };
    match bcftools::stats::stats_native(&input) {
        Ok(s) => {
            println!("# rubam-bcftools stats");
            println!("# input: {input}");
            println!("total_records\t{}", s.total_records);
            println!("snps\t{}", s.snps);
            println!("indels\t{}", s.indels);
            println!("mnps\t{}", s.mnps);
            println!("complex\t{}", s.complex);
            println!("transitions\t{}", s.transitions);
            println!("transversions\t{}", s.transversions);
            let ratio = if s.transversions == 0 {
                0.0
            } else {
                s.transitions as f64 / s.transversions as f64
            };
            println!("ts_tv_ratio\t{:.4}", ratio);
            println!("# per-sample (sample, hom_ref, het, hom_alt, missing)");
            for (name, ss) in &s.samples {
                println!(
                    "sample\t{name}\t{}\t{}\t{}\t{}",
                    ss.hom_ref, ss.het, ss.hom_alt, ss.missing
                );
            }
            0
        }
        Err(e) => {
            eprintln!("stats: {e}");
            1
        }
    }
}
