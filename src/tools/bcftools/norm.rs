//! `bcftools norm` — split multi-allelic + left-align indels with FASTA.
//!
//! v0.3 MVP implements:
//!   - `-m -`: split multi-allelic sites into one ALT per record.
//!   - `-f/--reference ref.fa`: left-align indels (simple last-base-match
//!     algorithm; MNV/complex pass through unchanged).
//!
//! `-m +` (join) and `-a` (atomize) land in v0.3.x.

use std::fs;
use std::io::{self, BufRead, BufWriter};
use std::path::Path;

use noodles::bcf;
use noodles::bgzf;
use noodles::core::Position;
use noodles::fasta;
use noodles::vcf;
use noodles::vcf::variant::record_buf::samples::sample::value::{
    genotype::Allele, Array, Genotype, Value as SampleValue,
};
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
// Writer enum (same pattern as concat.rs)
// ---------------------------------------------------------------------------

enum Writer {
    Plain(vcf::io::Writer<BufWriter<fs::File>>),
    Bgzf(vcf::io::Writer<bgzf::io::Writer<BufWriter<fs::File>>>),
    Bcf(bcf::io::Writer<bgzf::io::Writer<BufWriter<fs::File>>>),
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
// Input reader
// ---------------------------------------------------------------------------

fn read_all_records(input: &str) -> io::Result<(vcf::Header, Vec<vcf::variant::RecordBuf>)> {
    if input.ends_with(".bcf") {
        let f = fs::File::open(input)?;
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
// Multi-allelic split helpers
// ---------------------------------------------------------------------------

/// Remap a single allele index for split index `alt_i` (0-based among ALTs).
/// - 0 (ref) → 0
/// - None (missing) → None
/// - alt_i + 1 → 1
/// - anything else → 0  (degraded per spec comment in plan)
fn remap_allele_index(idx: Option<usize>, alt_i: usize) -> Option<usize> {
    match idx {
        None => None,
        Some(0) => Some(0),
        Some(n) if n == alt_i + 1 => Some(1),
        Some(_) => Some(0),
    }
}

/// Slice a FORMAT `Number=A` array value to keep only index `alt_i`.
fn slice_array_a(val: &SampleValue, alt_i: usize) -> SampleValue {
    match val {
        SampleValue::Array(Array::Integer(v)) => SampleValue::from(vec![v
            .get(alt_i)
            .copied()
            .flatten()
            .map(Some)
            .unwrap_or(None)]),
        SampleValue::Array(Array::Float(v)) => SampleValue::from(vec![v
            .get(alt_i)
            .copied()
            .flatten()
            .map(Some)
            .unwrap_or(None)]),
        SampleValue::Array(Array::Character(v)) => SampleValue::from(vec![v
            .get(alt_i)
            .copied()
            .flatten()
            .map(Some)
            .unwrap_or(None)]),
        SampleValue::Array(Array::String(v)) => SampleValue::from(vec![v
            .get(alt_i)
            .cloned()
            .flatten()
            .map(Some)
            .unwrap_or(None)]),
        other => other.clone(),
    }
}

/// Slice a FORMAT `Number=R` array value to keep ref (index 0) + index `alt_i+1`.
fn slice_array_r(val: &SampleValue, alt_i: usize) -> SampleValue {
    match val {
        SampleValue::Array(Array::Integer(v)) => {
            let ref_val = v.first().copied().flatten();
            let alt_val = v.get(alt_i + 1).copied().flatten();
            SampleValue::from(vec![
                ref_val.map(Some).unwrap_or(None),
                alt_val.map(Some).unwrap_or(None),
            ])
        }
        SampleValue::Array(Array::Float(v)) => {
            let ref_val = v.first().copied().flatten();
            let alt_val = v.get(alt_i + 1).copied().flatten();
            SampleValue::from(vec![
                ref_val.map(Some).unwrap_or(None),
                alt_val.map(Some).unwrap_or(None),
            ])
        }
        SampleValue::Array(Array::Character(v)) => {
            let ref_val = v.first().copied().flatten();
            let alt_val = v.get(alt_i + 1).copied().flatten();
            SampleValue::from(vec![
                ref_val.map(Some).unwrap_or(None),
                alt_val.map(Some).unwrap_or(None),
            ])
        }
        SampleValue::Array(Array::String(v)) => {
            let ref_val = v.first().cloned().flatten();
            let alt_val = v.get(alt_i + 1).cloned().flatten();
            SampleValue::from(vec![
                ref_val.map(Some).unwrap_or(None),
                alt_val.map(Some).unwrap_or(None),
            ])
        }
        other => other.clone(),
    }
}

/// Build one split record from `record` by keeping only ALT at 0-based index `alt_i`.
fn split_one(
    record: &vcf::variant::RecordBuf,
    header: &vcf::Header,
    alt_i: usize,
    alt_str: &str,
) -> vcf::variant::RecordBuf {
    use noodles::vcf::header::record::value::map::format::Number as FormatNumber;
    use noodles::vcf::header::record::value::map::info::Number as InfoNumber;
    use noodles::vcf::variant::record_buf::{AlternateBases, Samples};

    // --- Build new ALT bases ---
    let new_alts = AlternateBases::from(vec![alt_str.to_string()]);

    // --- Rebuild samples ---
    // Decompose via the From impl so we can access the raw value rows.
    let (keys, raw_values): (_, Vec<Vec<Option<SampleValue>>>) = record.samples().clone().into();
    let mut new_values: Vec<Vec<Option<SampleValue>>> = Vec::with_capacity(raw_values.len());

    for sample_row in &raw_values {
        let mut new_row: Vec<Option<SampleValue>> = Vec::with_capacity(sample_row.len());

        for (col_idx, val) in sample_row.iter().enumerate() {
            let key_name: Option<&str> = keys.as_ref().get_index(col_idx).map(|s| s.as_str());

            let new_val: Option<SampleValue> = match (key_name, val) {
                // GT: remap allele indices
                (Some("GT"), Some(SampleValue::Genotype(gt))) => {
                    let new_alleles: Vec<Allele> = gt
                        .as_ref()
                        .iter()
                        .map(|a| {
                            let new_pos = remap_allele_index(a.position(), alt_i);
                            Allele::new(new_pos, a.phasing())
                        })
                        .collect();
                    Some(SampleValue::Genotype(Genotype::from_iter(new_alleles)))
                }
                (Some(key), Some(inner)) => {
                    // Look up Number in the FORMAT header definitions.
                    let fmt_def = header.formats().get(key);
                    let number = fmt_def.map(|d| d.number());
                    match number {
                        Some(FormatNumber::AlternateBases) => Some(slice_array_a(inner, alt_i)),
                        Some(FormatNumber::ReferenceAlternateBases) => {
                            Some(slice_array_r(inner, alt_i))
                        }
                        _ => Some(inner.clone()),
                    }
                }
                (_, None) => None,
                (None, Some(inner)) => Some(inner.clone()),
            };
            new_row.push(new_val);
        }

        new_values.push(new_row);
    }

    let new_samples = Samples::new(keys, new_values);

    // --- Build new INFO: handle Number=A and Number=R ---
    let mut new_info = record.info().clone();
    for (k, v) in new_info.as_mut().iter_mut() {
        // Look up Number in header INFO definitions
        let info_def = header.infos().get(k.as_str());
        let number = info_def.map(|d| d.number());
        if let Some(val) = v {
            match number {
                Some(InfoNumber::AlternateBases) => {
                    // Keep only index alt_i
                    let sliced = slice_info_value_a(val, alt_i);
                    *val = sliced;
                }
                Some(InfoNumber::ReferenceAlternateBases) => {
                    let sliced = slice_info_value_r(val, alt_i);
                    *val = sliced;
                }
                _ => {}
            }
        }
    }

    let mut out = record.clone();
    *out.alternate_bases_mut() = new_alts;
    *out.samples_mut() = new_samples;
    *out.info_mut() = new_info;
    out
}

/// Slice an INFO `Number=A` value to keep only index `alt_i`.
fn slice_info_value_a(
    val: &noodles::vcf::variant::record_buf::info::field::Value,
    alt_i: usize,
) -> noodles::vcf::variant::record_buf::info::field::Value {
    use noodles::vcf::variant::record_buf::info::field::value::Array;
    use noodles::vcf::variant::record_buf::info::field::Value;
    match val {
        Value::Array(Array::Integer(v)) => Value::Array(Array::Integer(vec![v
            .get(alt_i)
            .copied()
            .flatten()
            .map(Some)
            .unwrap_or(None)])),
        Value::Array(Array::Float(v)) => Value::Array(Array::Float(vec![v
            .get(alt_i)
            .copied()
            .flatten()
            .map(Some)
            .unwrap_or(None)])),
        Value::Array(Array::Character(v)) => Value::Array(Array::Character(vec![v
            .get(alt_i)
            .copied()
            .flatten()
            .map(Some)
            .unwrap_or(None)])),
        Value::Array(Array::String(v)) => Value::Array(Array::String(vec![v
            .get(alt_i)
            .cloned()
            .flatten()
            .map(Some)
            .unwrap_or(None)])),
        other => other.clone(),
    }
}

/// Slice an INFO `Number=R` value to keep ref (index 0) + index `alt_i+1`.
fn slice_info_value_r(
    val: &noodles::vcf::variant::record_buf::info::field::Value,
    alt_i: usize,
) -> noodles::vcf::variant::record_buf::info::field::Value {
    use noodles::vcf::variant::record_buf::info::field::value::Array;
    use noodles::vcf::variant::record_buf::info::field::Value;
    match val {
        Value::Array(Array::Integer(v)) => {
            let r = v.first().copied().flatten();
            let a = v.get(alt_i + 1).copied().flatten();
            Value::Array(Array::Integer(vec![
                r.map(Some).unwrap_or(None),
                a.map(Some).unwrap_or(None),
            ]))
        }
        Value::Array(Array::Float(v)) => {
            let r = v.first().copied().flatten();
            let a = v.get(alt_i + 1).copied().flatten();
            Value::Array(Array::Float(vec![
                r.map(Some).unwrap_or(None),
                a.map(Some).unwrap_or(None),
            ]))
        }
        Value::Array(Array::Character(v)) => {
            let r = v.first().copied().flatten();
            let a = v.get(alt_i + 1).copied().flatten();
            Value::Array(Array::Character(vec![
                r.map(Some).unwrap_or(None),
                a.map(Some).unwrap_or(None),
            ]))
        }
        Value::Array(Array::String(v)) => {
            let r = v.first().cloned().flatten();
            let a = v.get(alt_i + 1).cloned().flatten();
            Value::Array(Array::String(vec![
                r.map(Some).unwrap_or(None),
                a.map(Some).unwrap_or(None),
            ]))
        }
        other => other.clone(),
    }
}

// ---------------------------------------------------------------------------
// Left-align indels
// ---------------------------------------------------------------------------

/// Open an indexed FASTA reader, building the .fai if it doesn't exist yet.
fn open_fasta(
    ref_path: &Path,
) -> io::Result<fasta::io::IndexedReader<fasta::io::BufReader<fs::File>>> {
    let fai_path = {
        let mut p = ref_path.as_os_str().to_owned();
        p.push(".fai");
        std::path::PathBuf::from(p)
    };

    if !fai_path.exists() {
        let idx = fasta::fs::index(ref_path)?;
        // Write the .fai next to the FASTA
        let mut fai_writer =
            fasta::fai::io::Writer::new(BufWriter::new(fs::File::create(&fai_path)?));
        fai_writer.write_index(&idx)?;
    }

    fasta::io::indexed_reader::Builder::default().build_from_path(ref_path)
}

/// Fetch a single base (1-based) from the FASTA. Returns None on any error.
fn fasta_base(
    reader: &mut fasta::io::IndexedReader<fasta::io::BufReader<fs::File>>,
    chrom: &str,
    pos1: usize, // 1-based
) -> Option<u8> {
    let start = Position::try_from(pos1).ok()?;
    let region = noodles::core::Region::new(chrom, start..=start);
    let record = reader.query(&region).ok()?;
    let seq = record.sequence();
    seq.as_ref().first().copied()
}

/// Try to left-align an indel. Returns (new_pos, new_ref, new_alt, shifted).
/// `pos` is 1-based.
fn left_align(
    reader: &mut fasta::io::IndexedReader<fasta::io::BufReader<fs::File>>,
    chrom: &str,
    pos: usize, // 1-based
    ref_bases: &str,
    alt_bases: &str,
) -> (usize, String, String, bool) {
    let mut p = pos;
    let mut r: Vec<u8> = ref_bases.as_bytes().to_vec();
    let mut a: Vec<u8> = alt_bases.as_bytes().to_vec();
    let mut shifted = false;

    // Simple algorithm: while last base of REF == last base of ALT AND p > 1,
    // fetch the base at p-1, shift left.
    loop {
        // Need at least 2 chars in both (the anchor base must remain)
        if r.is_empty() || a.is_empty() {
            break;
        }
        if *r.last().unwrap() != *a.last().unwrap() {
            break;
        }
        if p <= 1 {
            break;
        }
        let prior = match fasta_base(reader, chrom, p - 1) {
            Some(b) => b.to_ascii_uppercase(),
            None => break,
        };
        // Shift: prepend prior, drop last base
        r.pop();
        r.insert(0, prior);
        a.pop();
        a.insert(0, prior);
        p -= 1;
        shifted = true;
    }

    // Trim common leading prefix (keep at least 1 anchor base)
    while r.len() > 1 && a.len() > 1 && r[0] == a[0] {
        r.remove(0);
        a.remove(0);
        p += 1;
    }

    (
        p,
        String::from_utf8_lossy(&r).into_owned(),
        String::from_utf8_lossy(&a).into_owned(),
        shifted,
    )
}

// ---------------------------------------------------------------------------
// Core native function
// ---------------------------------------------------------------------------

/// Native norm. Returns (records_in, records_out, left_aligned_count).
pub fn norm_native(
    input: &str,
    output: &Path,
    format: OutputFormat,
    split_multiallelic: bool,
    reference: Option<&Path>,
) -> io::Result<(usize, usize, usize)> {
    // 1. Read all records.
    let (header, records) = read_all_records(input)?;

    // 2. Optionally open FASTA reader.
    let mut fasta_reader = reference.map(open_fasta).transpose()?;

    // 3. Open output writer.
    let f = fs::File::create(output)?;
    let buf = BufWriter::new(f);
    let mut writer: Writer = match format {
        OutputFormat::Vcf => Writer::Plain(vcf::io::Writer::new(buf)),
        OutputFormat::VcfGz => Writer::Bgzf(vcf::io::Writer::new(bgzf::io::Writer::new(buf))),
        OutputFormat::Bcf => Writer::Bcf(bcf::io::Writer::new(buf)),
    };
    writer.write_header(&header)?;

    let records_in = records.len();
    let mut records_out = 0usize;
    let mut left_aligned = 0usize;

    for record in &records {
        // Collect ALT strings.
        let alts: Vec<String> = record.alternate_bases().as_ref().iter().cloned().collect();

        // Decide whether to split.
        if split_multiallelic && alts.len() > 1 {
            for (i, alt_str) in alts.iter().enumerate() {
                let mut split_rec = split_one(record, &header, i, alt_str);

                // Optionally left-align the split record.
                if let Some(ref mut fr) = fasta_reader {
                    let chrom = split_rec.reference_sequence_name().to_string();
                    let pos = split_rec.variant_start().map(|p| p.get()).unwrap_or(1);
                    let ref_b = split_rec.reference_bases().to_string();
                    let alt_b = split_rec
                        .alternate_bases()
                        .as_ref()
                        .first()
                        .cloned()
                        .unwrap_or_default();

                    if ref_b.len() != alt_b.len() {
                        let (new_pos, new_ref, new_alt, shifted) =
                            left_align(fr, &chrom, pos, &ref_b, &alt_b);
                        if shifted || new_pos != pos {
                            *split_rec.variant_start_mut() = Position::try_from(new_pos).ok();
                            *split_rec.reference_bases_mut() = new_ref;
                            *split_rec.alternate_bases_mut() =
                                vcf::variant::record_buf::AlternateBases::from(vec![new_alt]);
                            left_aligned += 1;
                        }
                    }
                }

                writer.write_record(&header, &split_rec)?;
                records_out += 1;
            }
        } else {
            // Pass-through (possibly with left-align for single-ALT indels).
            let mut out_rec = record.clone();

            if let Some(ref mut fr) = fasta_reader {
                let chrom = out_rec.reference_sequence_name().to_string();
                let pos = out_rec.variant_start().map(|p| p.get()).unwrap_or(1);
                let ref_b = out_rec.reference_bases().to_string();
                let alts_single: Vec<String> =
                    out_rec.alternate_bases().as_ref().iter().cloned().collect();

                if alts_single.len() == 1 {
                    let alt_b = &alts_single[0];
                    if ref_b.len() != alt_b.len() {
                        let (new_pos, new_ref, new_alt, shifted) =
                            left_align(fr, &chrom, pos, &ref_b, alt_b);
                        if shifted || new_pos != pos {
                            *out_rec.variant_start_mut() = Position::try_from(new_pos).ok();
                            *out_rec.reference_bases_mut() = new_ref;
                            *out_rec.alternate_bases_mut() =
                                vcf::variant::record_buf::AlternateBases::from(vec![new_alt]);
                            left_aligned += 1;
                        }
                    }
                }
                // Multi-allelic with no split flag: don't try to left-align (skip).
            }

            writer.write_record(&header, &out_rec)?;
            records_out += 1;
        }
    }

    writer.finish()?;
    Ok((records_in, records_out, left_aligned))
}

// ---------------------------------------------------------------------------
// Python-facing wrapper
// ---------------------------------------------------------------------------

#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(name = "bcftools_norm")]
#[pyo3(signature = (
    input,
    output,
    *,
    output_type = "v",
    multiallelic = "",
    reference = None,
))]
pub fn bcftools_norm(
    input: &str,
    output: &str,
    output_type: &str,
    multiallelic: &str,
    reference: Option<&str>,
) -> PyResult<PyObject> {
    let format = OutputFormat::from_str(output_type)
        .ok_or_else(|| PyValueError::new_err(format!("unknown output_type {output_type:?}")))?;
    let split = match multiallelic {
        "" => false,
        "-" => true,
        "+" => {
            return Err(PyValueError::new_err(
                "-m + (join) lands in v0.3.x; use -m -",
            ))
        }
        other => return Err(PyValueError::new_err(format!("unknown -m {other:?}"))),
    };
    let ref_path = reference.map(Path::new);
    let (in_n, out_n, la_n) = norm_native(input, Path::new(output), format, split, ref_path)
        .map_err(|e| PyIOError::new_err(format!("bcftools norm: {e}")))?;
    Python::with_gil(|py| {
        let dict = pyo3::types::PyDict::new(py);
        dict.set_item("records_in", in_n)?;
        dict.set_item("records_out", out_n)?;
        dict.set_item("left_aligned", la_n)?;
        Ok(dict.into())
    })
}
