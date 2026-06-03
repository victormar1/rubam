//! `bcftools stats` — per-sample summary (Ts/Tv, het/hom, totals).

use std::io::{self, BufRead};

use noodles::bcf;
use noodles::vcf;
#[cfg(feature = "python")]
use pyo3::exceptions::PyIOError;
#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::types::PyDict;

// ---------------------------------------------------------------------------
// Result types
// ---------------------------------------------------------------------------

#[derive(Default, Debug)]
pub struct StatsResult {
    pub total_records: usize,
    pub snps: usize,
    pub indels: usize,
    pub mnps: usize,
    pub complex: usize,
    pub transitions: usize,
    pub transversions: usize,
    pub samples: Vec<(String, SampleStats)>,
}

#[derive(Default, Debug, Clone)]
pub struct SampleStats {
    pub hom_ref: usize,
    pub het: usize,
    pub hom_alt: usize,
    pub missing: usize,
}

// ---------------------------------------------------------------------------
// Variant classification
// ---------------------------------------------------------------------------

#[derive(Debug, PartialEq, Eq)]
enum Class {
    Snp,
    Indel,
    Mnp,
    Complex,
}

fn classify_record(ref_bases: &str, alt_bases: &[String]) -> Class {
    if alt_bases.is_empty() {
        return Class::Complex;
    }
    let ref_len = ref_bases.len();
    let all_one = ref_len == 1 && alt_bases.iter().all(|a| a.len() == 1);
    if all_one {
        return Class::Snp;
    }
    let any_indel = alt_bases.iter().any(|a| a.len() != ref_len);
    if any_indel {
        return Class::Indel;
    }
    let all_same_len = alt_bases.iter().all(|a| a.len() == ref_len);
    if ref_len > 1 && all_same_len {
        return Class::Mnp;
    }
    Class::Complex
}

fn is_transition(r: u8, a: u8) -> bool {
    matches!(
        (r, a),
        (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C')
    )
}

// ---------------------------------------------------------------------------
// Sample GT counting
// ---------------------------------------------------------------------------

/// Classify a GT from an iterator of `(Option<usize>, Phasing)` results.
/// Returns (hom_ref_delta, het_delta, hom_alt_delta, missing_delta).
fn classify_gt<I>(mut allele_iter: I) -> (usize, usize, usize, usize)
where
    I: Iterator<
        Item = io::Result<(
            Option<usize>,
            vcf::variant::record::samples::series::value::genotype::Phasing,
        )>,
    >,
{
    let mut indices: Vec<Option<usize>> = Vec::new();
    for result in &mut allele_iter {
        match result {
            Ok((idx, _phasing)) => indices.push(idx),
            Err(_) => {
                // Malformed genotype — treat as missing.
                return (0, 0, 0, 1);
            }
        }
    }

    if indices.is_empty() {
        return (0, 0, 0, 1);
    }

    // If any allele is None → missing
    if indices.iter().any(|i| i.is_none()) {
        return (0, 0, 0, 1);
    }

    let vals: Vec<usize> = indices.into_iter().map(|i| i.unwrap()).collect();

    // All zero → hom_ref
    if vals.iter().all(|&v| v == 0) {
        return (1, 0, 0, 0);
    }
    // All same non-zero → hom_alt
    if vals.iter().all(|&v| v == vals[0]) && vals[0] != 0 {
        return (0, 0, 1, 0);
    }
    // Otherwise → het
    (0, 1, 0, 0)
}

// ---------------------------------------------------------------------------
// Input reader
// ---------------------------------------------------------------------------

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
// Native stats
// ---------------------------------------------------------------------------

/// Native stats. Streams the input.
pub fn stats_native(input: &str) -> io::Result<StatsResult> {
    use noodles::vcf::variant::record::samples::Sample as _;
    use noodles::vcf::variant::record::AlternateBases as _;

    let (header, records) = read_all_records(input)?;

    // Initialize per-sample counters from header.
    let sample_names: Vec<String> = header.sample_names().iter().cloned().collect();

    let mut result = StatsResult {
        samples: sample_names
            .iter()
            .map(|n| (n.clone(), SampleStats::default()))
            .collect(),
        ..Default::default()
    };

    for record in &records {
        result.total_records += 1;

        // Collect ALT strings.
        let ref_bases = record.reference_bases();
        let alt_bases: Vec<String> = record
            .alternate_bases()
            .iter()
            .filter_map(|r| r.ok().map(|s| s.to_string()))
            .collect();

        // Classify the record.
        let class = classify_record(ref_bases, &alt_bases);
        match class {
            Class::Snp => {
                result.snps += 1;
                // Count Ts/Tv for each ALT allele.
                let ref_upper = ref_bases.to_ascii_uppercase();
                if ref_upper.len() == 1 {
                    let r = ref_upper.as_bytes()[0];
                    for alt in &alt_bases {
                        let alt_upper = alt.to_ascii_uppercase();
                        if alt_upper.len() == 1 {
                            let a = alt_upper.as_bytes()[0];
                            if is_transition(r, a) {
                                result.transitions += 1;
                            } else {
                                result.transversions += 1;
                            }
                        }
                    }
                }
            }
            Class::Indel => result.indels += 1,
            Class::Mnp => result.mnps += 1,
            Class::Complex => result.complex += 1,
        }

        // Per-sample GT counting.
        let samples_data = record.samples();
        for (sample_idx, sample_name) in sample_names.iter().enumerate() {
            let sample_opt = samples_data.get(&header, sample_name);
            let stat = &mut result.samples[sample_idx].1;

            let (hr, het, ha, mis) = match sample_opt {
                None => (0, 0, 0, 1),
                Some(sample) => {
                    // Find GT field by iterating FORMAT fields.
                    let mut gt_result = (0usize, 0usize, 0usize, 1usize); // default: missing
                    for field_result in sample.iter(&header) {
                        match field_result {
                            Ok(("GT", Some(v))) => {
                                use vcf::variant::record::samples::series::Value;
                                if let Value::Genotype(gt) = v {
                                    gt_result = classify_gt(gt.iter());
                                }
                                break;
                            }
                            Ok(("GT", None)) => {
                                // GT present but missing value.
                                gt_result = (0, 0, 0, 1);
                                break;
                            }
                            _ => {}
                        }
                    }
                    gt_result
                }
            };

            stat.hom_ref += hr;
            stat.het += het;
            stat.hom_alt += ha;
            stat.missing += mis;
        }
    }

    Ok(result)
}

// ---------------------------------------------------------------------------
// Python-facing wrapper
// ---------------------------------------------------------------------------

#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(name = "bcftools_stats")]
pub fn bcftools_stats(py: Python<'_>, input: &str) -> PyResult<PyObject> {
    let s = stats_native(input).map_err(|e| PyIOError::new_err(format!("bcftools stats: {e}")))?;
    let dict = PyDict::new(py);
    dict.set_item("total_records", s.total_records)?;
    dict.set_item("snps", s.snps)?;
    dict.set_item("indels", s.indels)?;
    dict.set_item("mnps", s.mnps)?;
    dict.set_item("complex", s.complex)?;
    dict.set_item("transitions", s.transitions)?;
    dict.set_item("transversions", s.transversions)?;
    dict.set_item(
        "ts_tv_ratio",
        if s.transversions == 0 {
            0.0_f64
        } else {
            s.transitions as f64 / s.transversions as f64
        },
    )?;
    let samples = PyDict::new(py);
    for (name, ss) in &s.samples {
        let sd = PyDict::new(py);
        sd.set_item("hom_ref", ss.hom_ref)?;
        sd.set_item("het", ss.het)?;
        sd.set_item("hom_alt", ss.hom_alt)?;
        sd.set_item("missing", ss.missing)?;
        samples.set_item(name, sd)?;
    }
    dict.set_item("samples", samples)?;
    Ok(dict.into())
}
