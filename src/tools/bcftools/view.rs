//! `bcftools view` — region/sample filter + output.

use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufWriter};
use std::path::Path;

use noodles::bcf;
use noodles::bgzf;
use noodles::vcf;
use noodles::vcf::variant::record_buf::samples::Keys as SampleKeys;
#[cfg(feature = "python")]
use pyo3::exceptions::{PyIOError, PyValueError};
#[cfg(feature = "python")]
use pyo3::prelude::*;

// ---------------------------------------------------------------------------
// Output format
// ---------------------------------------------------------------------------

/// Output format selector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    Vcf,
    VcfGz,
    Bcf,
}

impl OutputFormat {
    pub fn from_str(s: &str) -> Option<Self> {
        match s {
            "v" | "vcf" => Some(Self::Vcf),
            "z" | "vcf.gz" => Some(Self::VcfGz),
            "b" | "bcf" => Some(Self::Bcf),
            _ => None,
        }
    }
}

// ---------------------------------------------------------------------------
// Region parsing
// ---------------------------------------------------------------------------

/// Parse `"chr:start-end"` (1-based inclusive) or `"chr"` (full contig).
fn parse_region(s: &str) -> io::Result<(String, Option<(usize, usize)>)> {
    if let Some((chrom, range)) = s.split_once(':') {
        let (start_str, end_str) = range.split_once('-').ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "region must be chr:start-end")
        })?;
        let start: usize = start_str.parse().map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidInput, "region start not a number")
        })?;
        let end: usize = end_str
            .parse()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "region end not a number"))?;
        Ok((chrom.to_string(), Some((start, end))))
    } else {
        Ok((s.to_string(), None))
    }
}

/// Returns `true` if the record overlaps the region.
fn record_in_region(
    rec: &vcf::variant::RecordBuf,
    chrom: &str,
    range: Option<(usize, usize)>,
) -> bool {
    if rec.reference_sequence_name() != chrom {
        return false;
    }
    if let Some((start, end)) = range {
        // VCF POS is 1-based.
        match rec.variant_start() {
            Some(pos) => {
                let p = pos.get();
                // The variant overlaps [start, end] when p <= end and p >= start.
                p >= start && p <= end
            }
            None => false,
        }
    } else {
        true
    }
}

// ---------------------------------------------------------------------------
// Sample-subset helpers
// ---------------------------------------------------------------------------

/// Build the list of sample indices to keep (in original order) given a keep-set.
fn keep_indices(header: &vcf::Header, keep: &HashSet<String>) -> Vec<usize> {
    header
        .sample_names()
        .iter()
        .enumerate()
        .filter_map(|(i, name)| {
            if keep.contains(name.as_str()) {
                Some(i)
            } else {
                None
            }
        })
        .collect()
}

/// Restrict the header to only the kept sample names (modifies in-place).
fn subset_header(header: &mut vcf::Header, keep: &HashSet<String>) {
    header.sample_names_mut().retain(|name| keep.contains(name));
}

/// Rebuild a `RecordBuf` Samples block keeping only the indicated sample indices.
fn subset_record_samples(record: &mut vcf::variant::RecordBuf, indices: &[usize]) {
    use vcf::variant::record_buf::samples::sample::Value;

    let orig = record.samples_mut();
    // Clone keys (FORMAT fields don't change, only which samples we keep).
    let new_keys: SampleKeys = orig.keys().clone();
    // Collect all per-sample value rows via the public iterator.
    // Sample::values() returns &[Option<Value>].
    let orig_rows: Vec<Vec<Option<Value>>> = orig.values().map(|s| s.values().to_vec()).collect();
    let new_values: Vec<Vec<Option<Value>>> = indices
        .iter()
        .filter_map(|&i| orig_rows.get(i).cloned())
        .collect();
    *orig = vcf::variant::record_buf::Samples::new(new_keys, new_values);
}

// ---------------------------------------------------------------------------
// Writer enum — avoids boxing and keeps try_finish() accessible
// ---------------------------------------------------------------------------

enum Writer {
    Plain(vcf::io::Writer<BufWriter<File>>),
    Bgzf(vcf::io::Writer<bgzf::io::Writer<BufWriter<File>>>),
    Bcf(bcf::io::Writer<bgzf::io::Writer<BufWriter<File>>>),
}

impl Writer {
    fn write_header(&mut self, header: &vcf::Header) -> io::Result<()> {
        match self {
            Self::Plain(w) => w.write_header(header),
            Self::Bgzf(w) => w.write_header(header),
            Self::Bcf(w) => w.write_header(header),
        }
    }

    fn write_record(
        &mut self,
        header: &vcf::Header,
        record: &vcf::variant::RecordBuf,
    ) -> io::Result<()> {
        use noodles::vcf::variant::io::Write as _;
        match self {
            Self::Plain(w) => w.write_variant_record(header, record),
            Self::Bgzf(w) => w.write_variant_record(header, record),
            Self::Bcf(w) => w.write_variant_record(header, record),
        }
    }

    fn finish(self) -> io::Result<()> {
        match self {
            Self::Plain(_) => Ok(()),
            Self::Bgzf(mut w) => w.get_mut().try_finish(),
            Self::Bcf(mut w) => w.try_finish(),
        }
    }
}

// ---------------------------------------------------------------------------
// Input — support VCF, VCF.gz and BCF via a unified RecordBuf iterator
// ---------------------------------------------------------------------------

fn read_all_records_vcf(path: &str) -> io::Result<(vcf::Header, Vec<vcf::variant::RecordBuf>)> {
    let mut reader = vcf::io::reader::Builder::default().build_from_path(path)?;
    let header = reader.read_header()?;
    let mut records: Vec<vcf::variant::RecordBuf> = Vec::new();
    for result in reader.record_bufs(&header) {
        records.push(result?);
    }
    Ok((header, records))
}

fn read_all_records_bcf(path: &str) -> io::Result<(vcf::Header, Vec<vcf::variant::RecordBuf>)> {
    let f = File::open(path)?;
    let mut reader = bcf::io::Reader::new(f);
    let header = reader.read_header()?;
    let mut buf = vcf::variant::RecordBuf::default();
    let mut records: Vec<vcf::variant::RecordBuf> = Vec::new();
    loop {
        match reader.read_record_buf(&header, &mut buf)? {
            0 => break,
            _ => records.push(buf.clone()),
        }
    }
    Ok((header, records))
}

fn read_input(input: &str) -> io::Result<(vcf::Header, Vec<vcf::variant::RecordBuf>)> {
    if input.ends_with(".bcf") {
        read_all_records_bcf(input)
    } else {
        read_all_records_vcf(input)
    }
}

// ---------------------------------------------------------------------------
// Core native function
// ---------------------------------------------------------------------------

/// Native (Rust-side) view. Returns the number of records written
/// (excluding header). `output` must be provided in v0.3 (stdout is deferred).
#[allow(clippy::too_many_arguments)]
pub fn view_native(
    input: &str,
    region: Option<&str>,
    samples: Option<&[String]>,
    output: Option<&Path>,
    format: OutputFormat,
    header_only: bool,
    no_header: bool,
) -> io::Result<usize> {
    // 1. Read everything into memory (RecordBuf).
    let (mut header, records) = read_input(input)?;

    // 2. Parse region filter.
    let region_filter: Option<(String, Option<(usize, usize)>)> =
        region.map(parse_region).transpose()?;

    // 3. Build sample-subset info (indices + modified header).
    let indices_opt: Option<Vec<usize>> = if let Some(keep_names) = samples {
        let keep_set: HashSet<String> = keep_names.iter().cloned().collect();
        let idx = keep_indices(&header, &keep_set);
        subset_header(&mut header, &keep_set);
        Some(idx)
    } else {
        None
    };

    // 4. Open output writer (requires header to be final before opening).
    let out_path = output.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "output path is required in v0.3",
        )
    })?;

    let mut writer: Writer = {
        let f = File::create(out_path)?;
        let buf = BufWriter::new(f);
        match format {
            OutputFormat::Vcf => Writer::Plain(vcf::io::Writer::new(buf)),
            OutputFormat::VcfGz => Writer::Bgzf(vcf::io::Writer::new(bgzf::io::Writer::new(buf))),
            OutputFormat::Bcf => Writer::Bcf(bcf::io::Writer::new(buf)),
        }
    };

    // 5. Write header unless --no-header.
    if !no_header {
        writer.write_header(&header)?;
    }

    // 6. Return early if --header-only.
    if header_only {
        writer.finish()?;
        return Ok(0);
    }

    // 7. Iterate records, apply filters, write.
    let mut count = 0usize;
    for mut rec in records {
        // Region filter.
        if let Some((ref chrom, range)) = region_filter {
            if !record_in_region(&rec, chrom, range) {
                continue;
            }
        }

        // Sample subset.
        if let Some(ref indices) = indices_opt {
            subset_record_samples(&mut rec, indices);
        }

        writer.write_record(&header, &rec)?;
        count += 1;
    }

    writer.finish()?;
    Ok(count)
}

// ---------------------------------------------------------------------------
// Python-facing wrapper
// ---------------------------------------------------------------------------

/// Python-facing wrapper.
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(name = "bcftools_view")]
#[pyo3(signature = (
    input,
    region = None,
    samples = None,
    output = None,
    output_type = "v",
    header_only = false,
    no_header = false,
))]
pub fn bcftools_view(
    input: &str,
    region: Option<&str>,
    samples: Option<Vec<String>>,
    output: Option<&str>,
    output_type: &str,
    header_only: bool,
    no_header: bool,
) -> PyResult<usize> {
    let format = OutputFormat::from_str(output_type)
        .ok_or_else(|| PyValueError::new_err(format!("unknown output_type {output_type:?}")))?;
    let out_path = output.map(Path::new);
    view_native(
        input,
        region,
        samples.as_deref(),
        out_path,
        format,
        header_only,
        no_header,
    )
    .map_err(|e| PyIOError::new_err(format!("bcftools view: {e}")))
}
