//! `bcftools query` — format-string extraction.

use std::fs::File;
use std::io::{self, BufRead, BufWriter, Write};
use std::path::Path;

use noodles::bcf;
use noodles::vcf;
#[cfg(feature = "python")]
use pyo3::exceptions::{PyIOError, PyValueError};
#[cfg(feature = "python")]
use pyo3::prelude::*;

// ---------------------------------------------------------------------------
// Token types
// ---------------------------------------------------------------------------

/// Top-level token in the format string.
#[derive(Debug)]
enum FormatToken {
    /// Raw literal bytes (escape sequences already resolved: `\t`, `\n`).
    Literal(String),
    Chrom,
    Pos,
    Id,
    Ref,
    Alt,
    Qual,
    Filter,
    /// `%INFO/KEY` → `Info("KEY")`.
    Info(String),
    /// `[...]` block — repeated once per sample.
    SampleLoop(Vec<SampleSubToken>),
}

/// Sub-token inside a `[...]` sample loop.
#[derive(Debug)]
enum SampleSubToken {
    Literal(String),
    /// `%SAMPLE` → sample name.
    Sample,
    /// Any other `%KEY` inside a sample loop → FORMAT field lookup.
    Field(String),
}

// ---------------------------------------------------------------------------
// Format string parser
// ---------------------------------------------------------------------------

/// Resolve `\t` / `\n` escape sequences in a literal chunk (kept for potential future use).
#[allow(dead_code)]
fn resolve_escapes(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut chars = s.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '\\' {
            match chars.peek() {
                Some('t') => {
                    chars.next();
                    out.push('\t');
                }
                Some('n') => {
                    chars.next();
                    out.push('\n');
                }
                _ => out.push('\\'),
            }
        } else {
            out.push(c);
        }
    }
    out
}

/// Parse a format string into a `Vec<FormatToken>`.
///
/// Raises an `io::Error` when nested `[...]` blocks are detected or when the
/// format contains a `[` without a matching `]`.
fn parse_format(fmt: &str) -> io::Result<Vec<FormatToken>> {
    let mut tokens: Vec<FormatToken> = Vec::new();
    let mut chars = fmt.chars().peekable();
    let mut lit = String::new();

    while let Some(c) = chars.next() {
        match c {
            // Escape sequences in literal context
            '\\' => match chars.peek() {
                Some('t') => {
                    chars.next();
                    lit.push('\t');
                }
                Some('n') => {
                    chars.next();
                    lit.push('\n');
                }
                _ => lit.push('\\'),
            },

            '%' => {
                // Flush pending literal.
                if !lit.is_empty() {
                    tokens.push(FormatToken::Literal(lit.drain(..).collect()));
                }
                // Collect identifier (alphanumeric + underscore + '/').
                let mut ident = String::new();
                while let Some(&nc) = chars.peek() {
                    if nc.is_alphanumeric() || nc == '_' || nc == '/' {
                        ident.push(nc);
                        chars.next();
                    } else {
                        break;
                    }
                }
                let tok = match ident.as_str() {
                    "CHROM" => FormatToken::Chrom,
                    "POS" => FormatToken::Pos,
                    "ID" => FormatToken::Id,
                    "REF" => FormatToken::Ref,
                    "ALT" => FormatToken::Alt,
                    "QUAL" => FormatToken::Qual,
                    "FILTER" => FormatToken::Filter,
                    _ if ident.starts_with("INFO/") => {
                        FormatToken::Info(ident["INFO/".len()..].to_string())
                    }
                    other => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            format!("unknown top-level placeholder %{other}"),
                        ))
                    }
                };
                tokens.push(tok);
            }

            '[' => {
                // Flush pending literal.
                if !lit.is_empty() {
                    tokens.push(FormatToken::Literal(lit.drain(..).collect()));
                }
                // Collect everything until matching `]`, then parse as sub-tokens.
                let mut inner = String::new();
                for nc in chars.by_ref() {
                    if nc == ']' {
                        break;
                    }
                    if nc == '[' {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "nested [ ] blocks are not supported",
                        ));
                    }
                    inner.push(nc);
                }
                tokens.push(FormatToken::SampleLoop(parse_sample_loop(&inner)?));
            }

            ']' => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "unexpected ']' outside a sample loop",
                ));
            }

            other => lit.push(other),
        }
    }

    if !lit.is_empty() {
        tokens.push(FormatToken::Literal(lit));
    }

    Ok(tokens)
}

/// Parse the inner content of a `[...]` block into `SampleSubToken`s.
fn parse_sample_loop(inner: &str) -> io::Result<Vec<SampleSubToken>> {
    let mut sub: Vec<SampleSubToken> = Vec::new();
    let mut chars = inner.chars().peekable();
    let mut lit = String::new();

    while let Some(c) = chars.next() {
        match c {
            '\\' => match chars.peek() {
                Some('t') => {
                    chars.next();
                    lit.push('\t');
                }
                Some('n') => {
                    chars.next();
                    lit.push('\n');
                }
                _ => lit.push('\\'),
            },

            '%' => {
                if !lit.is_empty() {
                    sub.push(SampleSubToken::Literal(lit.drain(..).collect()));
                }
                let mut ident = String::new();
                while let Some(&nc) = chars.peek() {
                    if nc.is_alphanumeric() || nc == '_' {
                        ident.push(nc);
                        chars.next();
                    } else {
                        break;
                    }
                }
                let tok = if ident == "SAMPLE" {
                    SampleSubToken::Sample
                } else {
                    SampleSubToken::Field(ident)
                };
                sub.push(tok);
            }

            other => lit.push(other),
        }
    }

    if !lit.is_empty() {
        sub.push(SampleSubToken::Literal(lit));
    }

    Ok(sub)
}

// ---------------------------------------------------------------------------
// Value rendering helpers
// ---------------------------------------------------------------------------

/// Render a `record_buf::info::field::Value` as text for query output.
fn render_info_value(v: &vcf::variant::record_buf::info::field::Value) -> String {
    use vcf::variant::record_buf::info::field::value::Array;
    use vcf::variant::record_buf::info::field::Value;
    match v {
        Value::Integer(n) => n.to_string(),
        Value::Float(f) => format!("{f}"),
        Value::Flag => "1".to_string(),
        Value::Character(c) => c.to_string(),
        Value::String(s) => s.clone(),
        Value::Array(arr) => match arr {
            Array::Integer(vals) => vals
                .iter()
                .map(|o| o.map_or(".".to_string(), |n| n.to_string()))
                .collect::<Vec<_>>()
                .join(","),
            Array::Float(vals) => vals
                .iter()
                .map(|o| o.map_or(".".to_string(), |f| f.to_string()))
                .collect::<Vec<_>>()
                .join(","),
            Array::Character(vals) => vals
                .iter()
                .map(|o| o.map_or(".".to_string(), |c| c.to_string()))
                .collect::<Vec<_>>()
                .join(","),
            Array::String(vals) => vals
                .iter()
                .map(|o| o.as_deref().unwrap_or(".").to_string())
                .collect::<Vec<_>>()
                .join(","),
        },
    }
}

/// Render a `record::samples::series::Value` (live/borrowed) as text.
fn render_sample_value<'a>(v: vcf::variant::record::samples::series::Value<'a>) -> String {
    use vcf::variant::record::samples::series::value::Array;
    use vcf::variant::record::samples::series::Value;
    match v {
        Value::Integer(n) => n.to_string(),
        Value::Float(f) => f.to_string(),
        Value::Character(c) => c.to_string(),
        Value::String(s) => s.to_string(),
        Value::Genotype(gt) => {
            // Render as "0/1", "./.", "1|0", etc.
            let mut out = String::new();
            for (i, r) in gt.iter().enumerate() {
                if let Ok((allele_idx, phasing)) = r {
                    use vcf::variant::record::samples::series::value::genotype::Phasing;
                    let sep = if i == 0 {
                        ""
                    } else {
                        match phasing {
                            Phasing::Phased => "|",
                            Phasing::Unphased => "/",
                        }
                    };
                    let a = match allele_idx {
                        Some(idx) => idx.to_string(),
                        None => ".".to_string(),
                    };
                    out.push_str(sep);
                    out.push_str(&a);
                }
            }
            out
        }
        Value::Array(arr) => match arr {
            Array::Integer(vals) => vals
                .iter()
                .map(|r| r.ok().flatten().map_or(".".to_string(), |n| n.to_string()))
                .collect::<Vec<_>>()
                .join(","),
            Array::Float(vals) => vals
                .iter()
                .map(|r| r.ok().flatten().map_or(".".to_string(), |f| f.to_string()))
                .collect::<Vec<_>>()
                .join(","),
            Array::Character(vals) => vals
                .iter()
                .map(|r| r.ok().flatten().map_or(".".to_string(), |c| c.to_string()))
                .collect::<Vec<_>>()
                .join(","),
            Array::String(vals) => vals
                .iter()
                .map(|r| r.ok().flatten().map_or(".".to_string(), |s| s.to_string()))
                .collect::<Vec<_>>()
                .join(","),
        },
    }
}

// ---------------------------------------------------------------------------
// Record expansion
// ---------------------------------------------------------------------------

/// Expand one `RecordBuf` record to text using the pre-parsed token list,
/// and write to `out`.
fn expand_record<W: Write>(
    out: &mut W,
    tokens: &[FormatToken],
    record: &vcf::variant::RecordBuf,
    header: &vcf::Header,
) -> io::Result<()> {
    use noodles::vcf::variant::record::samples::Sample as _;
    use noodles::vcf::variant::record::AlternateBases as _;
    use noodles::vcf::variant::record::Filters as _;
    use noodles::vcf::variant::record::Ids as _;

    for tok in tokens {
        match tok {
            FormatToken::Literal(s) => out.write_all(s.as_bytes())?,

            FormatToken::Chrom => {
                out.write_all(record.reference_sequence_name().as_bytes())?;
            }

            FormatToken::Pos => {
                let pos = record.variant_start().map(|p| p.get()).unwrap_or(0);
                write!(out, "{pos}")?;
            }

            FormatToken::Id => {
                let ids: Vec<String> = record.ids().iter().map(|s| s.to_string()).collect();
                if ids.is_empty() {
                    out.write_all(b".")?;
                } else {
                    out.write_all(ids.join(";").as_bytes())?;
                }
            }

            FormatToken::Ref => {
                out.write_all(record.reference_bases().as_bytes())?;
            }

            FormatToken::Alt => {
                let mut alts: Vec<String> = Vec::new();
                for r in record.alternate_bases().iter() {
                    alts.push(
                        r.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?
                            .to_string(),
                    );
                }
                if alts.is_empty() {
                    out.write_all(b".")?;
                } else {
                    out.write_all(alts.join(",").as_bytes())?;
                }
            }

            FormatToken::Qual => match record.quality_score() {
                Some(q) => write!(out, "{q}")?,
                None => out.write_all(b".")?,
            },

            FormatToken::Filter => {
                let mut filters: Vec<String> = Vec::new();
                for r in record.filters().iter(header) {
                    let f =
                        r.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
                    filters.push(f.to_string());
                }
                if filters.is_empty() {
                    out.write_all(b".")?;
                } else {
                    out.write_all(filters.join(";").as_bytes())?;
                }
            }

            FormatToken::Info(key) => match record.info().get(key.as_str()) {
                Some(Some(v)) => out.write_all(render_info_value(v).as_bytes())?,
                _ => out.write_all(b".")?,
            },

            FormatToken::SampleLoop(sub_tokens) => {
                let sample_names: Vec<&str> =
                    header.sample_names().iter().map(String::as_str).collect();
                let samples_data = record.samples();

                for (sample_idx, sample_name) in sample_names.iter().enumerate() {
                    // Get the sample object for this index.
                    let sample_opt = samples_data.get(header, sample_name);

                    for sub in sub_tokens {
                        match sub {
                            SampleSubToken::Literal(s) => out.write_all(s.as_bytes())?,

                            SampleSubToken::Sample => {
                                out.write_all(sample_name.as_bytes())?;
                            }

                            SampleSubToken::Field(field_key) => {
                                let rendered = match &sample_opt {
                                    None => ".".to_string(),
                                    Some(sample) => {
                                        // Walk FORMAT fields for this sample.
                                        let mut found: Option<String> = None;
                                        for r in sample.iter(header) {
                                            match r {
                                                Ok((k, Some(v))) if k == field_key.as_str() => {
                                                    found = Some(render_sample_value(v));
                                                    break;
                                                }
                                                Ok((k, None)) if k == field_key.as_str() => {
                                                    found = Some(".".to_string());
                                                    break;
                                                }
                                                _ => {}
                                            }
                                        }
                                        found.unwrap_or_else(|| ".".to_string())
                                    }
                                };
                                let _ = sample_idx; // suppress unused warning
                                out.write_all(rendered.as_bytes())?;
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Input reader (mirrors sort.rs)
// ---------------------------------------------------------------------------

fn read_all_records(input: &str) -> io::Result<(vcf::Header, Vec<vcf::variant::RecordBuf>)> {
    if input.ends_with(".bcf") {
        let f = File::open(input)?;
        let mut reader = bcf::io::Reader::new(f);
        let header = reader.read_header()?;
        let mut buf = vcf::variant::RecordBuf::default();
        let mut records: Vec<vcf::variant::RecordBuf> = Vec::new();
        loop {
            match reader.read_record_buf(&header, &mut buf) {
                Ok(0) => break,
                Ok(_) => records.push(buf.clone()),
                Err(e) => return Err(e),
            }
        }
        Ok((header, records))
    } else {
        let mut reader: vcf::io::Reader<Box<dyn BufRead>> =
            vcf::io::reader::Builder::default().build_from_path(input)?;
        let header = reader.read_header()?;
        let mut buf = vcf::variant::RecordBuf::default();
        let mut records: Vec<vcf::variant::RecordBuf> = Vec::new();
        loop {
            match reader.read_record_buf(&header, &mut buf) {
                Ok(0) => break,
                Ok(_) => records.push(buf.clone()),
                Err(e) => return Err(e),
            }
        }
        Ok((header, records))
    }
}

// ---------------------------------------------------------------------------
// Native entry point
// ---------------------------------------------------------------------------

/// Native query. Returns the number of records emitted.
pub fn query_native(input: &str, format: &str, output: &Path) -> io::Result<usize> {
    // 1. Parse the format string into tokens once.
    let tokens = parse_format(format)?;

    // 2. Read all records.
    let (header, records) = read_all_records(input)?;

    // 3. Open output.
    let f = File::create(output)?;
    let mut out = BufWriter::new(f);

    // 4. Expand each record.
    for record in &records {
        expand_record(&mut out, &tokens, record, &header)?;
    }

    out.flush()?;
    Ok(records.len())
}

// ---------------------------------------------------------------------------
// Python-facing wrapper
// ---------------------------------------------------------------------------

#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(name = "bcftools_query")]
#[pyo3(signature = (input, format, output))]
pub fn bcftools_query(input: &str, format: &str, output: &str) -> PyResult<usize> {
    if format.is_empty() {
        return Err(PyValueError::new_err("--format is required"));
    }
    query_native(input, format, Path::new(output))
        .map_err(|e| PyIOError::new_err(format!("bcftools query: {e}")))
}
