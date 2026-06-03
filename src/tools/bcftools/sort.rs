//! `bcftools sort` — in-memory chrom+pos sort.

use std::io::{self, BufRead, Write};
use std::path::Path;

use noodles::bcf;
use noodles::bgzf;
use noodles::vcf;
#[cfg(feature = "python")]
use pyo3::exceptions::{PyIOError, PyValueError};
#[cfg(feature = "python")]
use pyo3::prelude::*;

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

/// Return the declaration-order index of `name` in the header's contig list.
fn header_contig_index(header: &vcf::Header, name: &str) -> Option<usize> {
    header.contigs().get_index_of(name)
}

/// Native sort. Returns the number of records written.
pub fn sort_native(input: &str, output: &Path, format: OutputFormat) -> io::Result<usize> {
    // 1. Read all records + header (auto-detects VCF/VCF.gz/BCF).
    let (header, records) = read_all_records(input)?;

    // 2. Sort by (contig declaration index, position).
    let mut records = records;
    records.sort_by_key(|r| {
        let chrom = r.reference_sequence_name().to_string();
        let idx = header_contig_index(&header, &chrom).unwrap_or(usize::MAX);
        let pos = r.variant_start().map(|p| p.get()).unwrap_or(0);
        (idx, pos)
    });

    // 3. Write to output.
    let n = records.len();
    write_records(output, format, &header, &records)?;
    Ok(n)
}

/// Read header + all records from any of VCF / VCF.gz / BCF.
fn read_all_records(input: &str) -> io::Result<(vcf::Header, Vec<vcf::variant::RecordBuf>)> {
    if input.ends_with(".bcf") {
        let f = std::fs::File::open(input)?;
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
        // VCF or VCF.gz — Builder auto-detects compression.
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

/// Write header + records to `output` in the requested format.
fn write_records(
    output: &Path,
    format: OutputFormat,
    header: &vcf::Header,
    records: &[vcf::variant::RecordBuf],
) -> io::Result<()> {
    match format {
        OutputFormat::Vcf => {
            let f = std::fs::File::create(output)?;
            let mut writer: vcf::io::Writer<Box<dyn Write>> =
                vcf::io::Writer::new(Box::new(f) as Box<dyn Write>);
            writer.write_header(header)?;
            for r in records {
                use noodles::vcf::variant::io::Write as _;
                writer.write_variant_record(header, r)?;
            }
        }
        OutputFormat::VcfGz => {
            let f = std::fs::File::create(output)?;
            let bgzf_inner = bgzf::io::Writer::new(Box::new(f) as Box<dyn Write>);
            let mut writer: vcf::io::Writer<bgzf::io::Writer<Box<dyn Write>>> =
                vcf::io::Writer::new(bgzf_inner);
            writer.write_header(header)?;
            for r in records {
                use noodles::vcf::variant::io::Write as _;
                writer.write_variant_record(header, r)?;
            }
        }
        OutputFormat::Bcf => {
            let f = std::fs::File::create(output)?;
            let mut writer: bcf::io::Writer<bgzf::io::Writer<Box<dyn Write>>> =
                bcf::io::Writer::new(Box::new(f) as Box<dyn Write>);
            writer.write_header(header)?;
            for r in records {
                use noodles::vcf::variant::io::Write as _;
                writer.write_variant_record(header, r)?;
            }
        }
    }
    Ok(())
}

#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(name = "bcftools_sort")]
#[pyo3(signature = (input, output, output_type = "v"))]
pub fn bcftools_sort(input: &str, output: &str, output_type: &str) -> PyResult<usize> {
    let format = OutputFormat::from_str(output_type)
        .ok_or_else(|| PyValueError::new_err(format!("unknown output_type {output_type:?}")))?;
    sort_native(input, Path::new(output), format)
        .map_err(|e| PyIOError::new_err(format!("bcftools sort: {e}")))
}
