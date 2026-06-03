//! `bcftools concat` — glue multiple sorted VCFs / BCFs.

use std::fs::File;
use std::io::{self, BufWriter};
use std::path::Path;

use noodles::bcf;
use noodles::bgzf;
use noodles::vcf;
#[cfg(feature = "python")]
use pyo3::exceptions::{PyIOError, PyValueError};
#[cfg(feature = "python")]
use pyo3::prelude::*;

// ---------------------------------------------------------------------------
// Output format
// ---------------------------------------------------------------------------

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
// Input helpers
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
// Header compatibility check
// ---------------------------------------------------------------------------

fn headers_compatible(a: &vcf::Header, b: &vcf::Header) -> Option<String> {
    if a.sample_names() != b.sample_names() {
        return Some("sample names differ".to_string());
    }
    let a_contigs: Vec<&str> = a.contigs().iter().map(|(n, _)| n.as_str()).collect();
    let b_contigs: Vec<&str> = b.contigs().iter().map(|(n, _)| n.as_str()).collect();
    if a_contigs != b_contigs {
        return Some("contig list differs".to_string());
    }
    None
}

// ---------------------------------------------------------------------------
// Core native function
// ---------------------------------------------------------------------------

/// Native concat. Returns the total number of records written.
pub fn concat_native(inputs: &[String], output: &Path, format: OutputFormat) -> io::Result<usize> {
    if inputs.len() < 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "concat requires at least 2 input files",
        ));
    }

    // 1. Read all inputs.
    let mut all_data: Vec<(vcf::Header, Vec<vcf::variant::RecordBuf>)> = Vec::new();
    for path in inputs {
        all_data.push(read_input(path)?);
    }

    // 2. Clone the first header now — used for compatibility checks, writing,
    //    and record serialisation. Cloning here drops the borrow on all_data[0]
    //    before the subsequent iteration.
    let write_header = all_data[0].0.clone();

    for (i, (header, _)) in all_data.iter().enumerate().skip(1) {
        if let Some(reason) = headers_compatible(&write_header, header) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("input {i} header incompatible with first: {reason}"),
            ));
        }
    }

    // 3. Open output writer.
    let f = File::create(output)?;
    let buf = BufWriter::new(f);
    let mut writer: Writer = match format {
        OutputFormat::Vcf => Writer::Plain(vcf::io::Writer::new(buf)),
        OutputFormat::VcfGz => Writer::Bgzf(vcf::io::Writer::new(bgzf::io::Writer::new(buf))),
        OutputFormat::Bcf => Writer::Bcf(bcf::io::Writer::new(buf)),
    };

    // 4. Write header once.
    writer.write_header(&write_header)?;

    // 5. Stream records from all inputs.
    let mut count = 0usize;
    for (_, records) in &all_data {
        for rec in records {
            writer.write_record(&write_header, rec)?;
            count += 1;
        }
    }

    writer.finish()?;
    Ok(count)
}

// ---------------------------------------------------------------------------
// Python-facing wrapper
// ---------------------------------------------------------------------------

#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(name = "bcftools_concat")]
#[pyo3(signature = (inputs, output, output_type = "v"))]
pub fn bcftools_concat(inputs: Vec<String>, output: &str, output_type: &str) -> PyResult<usize> {
    let format = OutputFormat::from_str(output_type)
        .ok_or_else(|| PyValueError::new_err(format!("unknown output_type {output_type:?}")))?;
    concat_native(&inputs, Path::new(output), format)
        .map_err(|e| PyIOError::new_err(format!("bcftools concat: {e}")))
}
