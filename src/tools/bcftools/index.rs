//! `bcftools index` — TBI for VCF.gz, CSI for BCF / VCF.gz.

use std::fs::File;
use std::io::{self, BufReader};
use std::path::{Path, PathBuf};

use std::collections::HashMap;

use noodles::bcf;
use noodles::bgzf;
use noodles::csi;
use noodles::csi::binning_index::index::reference_sequence::index::BinnedIndex;
use noodles::tabix;
use noodles::vcf;
use noodles::vcf::variant::Record as _;
#[cfg(feature = "python")]
use pyo3::exceptions::{PyIOError, PyValueError};
#[cfg(feature = "python")]
use pyo3::prelude::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexKind {
    Tbi,
    Csi,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum InputKind {
    VcfGz,
    Bcf,
}

fn detect_format(path: &Path) -> io::Result<InputKind> {
    // BCF by extension
    if let Some(ext) = path.extension() {
        if ext == "bcf" {
            return Ok(InputKind::Bcf);
        }
    }
    // VCF.gz: name ends in ".vcf.gz"
    if path.to_string_lossy().ends_with(".vcf.gz") {
        return Ok(InputKind::VcfGz);
    }
    Err(io::Error::new(
        io::ErrorKind::InvalidInput,
        "input must be .vcf.gz or .bcf for indexing",
    ))
}

fn index_path(input: &Path, kind: IndexKind) -> PathBuf {
    let mut s = input.as_os_str().to_owned();
    match kind {
        IndexKind::Tbi => s.push(".tbi"),
        IndexKind::Csi => s.push(".csi"),
    }
    PathBuf::from(s)
}

/// Build a TBI index for a bgzipped VCF file.
fn build_tbi(input: &Path) -> io::Result<tabix::Index> {
    vcf::fs::index(input)
}

/// Build a CSI index for a bgzipped VCF file.
fn build_csi_vcfgz(input: &Path) -> io::Result<csi::Index> {
    use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
    use noodles::csi::binning_index::Indexer;

    let file = File::open(input).map(BufReader::new)?;
    let mut reader = vcf::io::Reader::new(bgzf::io::Reader::new(file));
    let header = reader.read_header()?;

    // Build a name→id map from the header contigs for CSI reference_sequence_id
    let contig_names: HashMap<String, usize> = header
        .contigs()
        .iter()
        .enumerate()
        .map(|(i, (name, _))| (name.to_string(), i))
        .collect();
    let contig_count = contig_names.len();

    let mut indexer = Indexer::<BinnedIndex>::default();
    let mut record = vcf::Record::default();
    let mut start_position = reader.get_ref().virtual_position();

    while reader.read_record(&mut record)? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        let ref_name = record.reference_sequence_name();
        let ref_id = contig_names.get(ref_name).copied().unwrap_or(contig_count);

        let start = record
            .variant_start()
            .transpose()?
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing variant start"))?;

        let end = record.variant_end(&header)?;

        indexer.add_record(Some((ref_id, start, end, true)), chunk)?;
        start_position = end_position;
    }

    Ok(indexer.build(header.contigs().len()))
}

/// Build a CSI index for a BCF file.
fn build_csi_bcf(input: &Path) -> io::Result<csi::Index> {
    bcf::fs::index(input)
}

/// Native indexer. Returns the path of the written index file.
pub fn index_native(input: &str, kind: IndexKind, force: bool) -> io::Result<PathBuf> {
    let path = Path::new(input);
    let fmt = detect_format(path)?;

    // Reject TBI on BCF
    if fmt == InputKind::Bcf && kind == IndexKind::Tbi {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "TBI index is not applicable to BCF files; use --csi",
        ));
    }

    let out_path = index_path(path, kind);

    // Guard: file already exists and force == false
    if out_path.exists() && !force {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!(
                "index file already exists: {}; use --force to overwrite",
                out_path.display()
            ),
        ));
    }

    match (fmt, kind) {
        (InputKind::VcfGz, IndexKind::Tbi) => {
            let idx = build_tbi(path)?;
            tabix::fs::write(&out_path, &idx)?;
        }
        (InputKind::VcfGz, IndexKind::Csi) => {
            let idx = build_csi_vcfgz(path)?;
            csi::fs::write(&out_path, &idx)?;
        }
        (InputKind::Bcf, IndexKind::Csi) => {
            let idx = build_csi_bcf(path)?;
            csi::fs::write(&out_path, &idx)?;
        }
        (InputKind::Bcf, IndexKind::Tbi) => {
            // Already rejected above, but keep exhaustive match
            unreachable!("BCF + TBI already rejected")
        }
    }

    Ok(out_path)
}

#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(name = "bcftools_index")]
#[pyo3(signature = (input, csi = false, force = false))]
pub fn bcftools_index(input: &str, csi: bool, force: bool) -> PyResult<String> {
    let kind = if csi { IndexKind::Csi } else { IndexKind::Tbi };
    let written = index_native(input, kind, force).map_err(|e| match e.kind() {
        io::ErrorKind::InvalidInput | io::ErrorKind::AlreadyExists => {
            PyValueError::new_err(format!("bcftools index: {e}"))
        }
        _ => PyIOError::new_err(format!("bcftools index: {e}")),
    })?;
    Ok(written.to_string_lossy().into_owned())
}
