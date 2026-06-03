//! rubam — Rust BAM library with two public surfaces.
//!
//!   - [`api`] — pure-Rust types ([`api::AlignmentFile`], [`api::AlignedSegment`],
//!     [`api::Header`], [`api::Cigar`], [`api::Aux`], [`api::Error`]). External
//!     Rust crates (e.g. HARMOS) depend on these directly without pulling in
//!     pyo3. Stable from v0.2.1.
//!
//!   - The `_rubam` Python extension module — pyo3-built classes for Python
//!     (`rubam.AlignmentFile`, `rubam.AlignedSegment`, …) plus the original
//!     batch helpers `get_depths`, `pileup_bases`, `count_reads`,
//!     `flag_stats`. Stable from v0.2.0.
//!
//! Both surfaces produce identical output on the same input, validated by the
//! existing pytest suite + the `tests/harmos_compat.rs` integration tests.
//! v0.2.2 will refactor the pyo3 layer to delegate to [`api`] so the noodles
//! call chain lives in one place; v0.2.1 keeps them as parallel implementations
//! to minimize the change-set for the HARMOS integration window.
//!
//! All BAM I/O goes through `noodles` — no `htslib`, no C dependency, builds
//! natively on `x86_64-pc-windows-msvc` without WSL or vcpkg.

#[cfg(feature = "python")]
use pyo3::prelude::*;

#[cfg(feature = "python")]
pub mod alignment;
pub mod api;
#[cfg(feature = "python")]
pub mod bgzf_file;
pub mod common;
pub mod depth;
#[cfg(feature = "python")]
pub mod fasta_file;
#[cfg(feature = "python")]
pub mod fastx_file;
#[cfg(feature = "python")]
pub mod pileup;
#[cfg(feature = "python")]
pub mod pileup_iter;
#[cfg(feature = "python")]
pub mod stats;
#[cfg(feature = "python")]
pub mod tabix_file;
pub mod tools;
#[cfg(feature = "python")]
pub mod variant;

/// Python module: `rubam._rubam`.
#[cfg(feature = "python")]
#[pymodule]
fn _rubam(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(depth::get_depths, m)?)?;
    m.add_function(wrap_pyfunction!(depth::get_depths_numpy, m)?)?;
    m.add_function(wrap_pyfunction!(depth::get_depths_regions, m)?)?;
    m.add_function(wrap_pyfunction!(pileup::pileup_bases, m)?)?;
    m.add_function(wrap_pyfunction!(stats::count_reads, m)?)?;
    m.add_function(wrap_pyfunction!(stats::flag_stats, m)?)?;
    m.add_function(wrap_pyfunction!(tools::sort::sort, m)?)?;
    m.add_function(wrap_pyfunction!(tools::index::index, m)?)?;
    m.add_function(wrap_pyfunction!(tools::view::view, m)?)?;
    m.add_function(wrap_pyfunction!(tools::merge::merge, m)?)?;
    m.add_function(wrap_pyfunction!(tools::idxstats::idxstats, m)?)?;
    m.add_function(wrap_pyfunction!(tools::flagstat::flagstat_py, m)?)?;
    m.add_function(wrap_pyfunction!(tools::faidx::faidx, m)?)?;
    m.add_function(wrap_pyfunction!(tools::calmd::calmd, m)?)?;
    m.add_function(wrap_pyfunction!(tools::bcftools::view::bcftools_view, m)?)?;
    m.add_function(wrap_pyfunction!(tools::bcftools::sort::bcftools_sort, m)?)?;
    m.add_function(wrap_pyfunction!(tools::bcftools::index::bcftools_index, m)?)?;
    m.add_function(wrap_pyfunction!(
        tools::bcftools::concat::bcftools_concat,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(tools::bcftools::query::bcftools_query, m)?)?;
    m.add_function(wrap_pyfunction!(tools::bcftools::stats::bcftools_stats, m)?)?;
    m.add_function(wrap_pyfunction!(tools::bcftools::norm::bcftools_norm, m)?)?;
    m.add_class::<alignment::AlignmentFile>()?;
    m.add_class::<alignment::AlignedSegment>()?;
    m.add_class::<alignment::AlignmentFileIter>()?;
    m.add_class::<alignment::AlignmentFileStreamIter>()?;
    m.add_class::<alignment::AlignmentFileFetchIter>()?;
    m.add_class::<alignment::Header>()?;
    m.add_class::<fasta_file::FastaFile>()?;
    m.add_class::<bgzf_file::BGZFile>()?;
    m.add_class::<fastx_file::FastxFile>()?;
    m.add_class::<fastx_file::FastxRecord>()?;
    m.add_class::<tabix_file::TabixFile>()?;
    m.add_class::<tabix_file::TabixFileIter>()?;
    m.add_class::<pileup_iter::PileupColumn>()?;
    m.add_class::<pileup_iter::PileupIter>()?;
    m.add_class::<variant::VariantFile>()?;
    m.add_class::<variant::VariantFileIter>()?;
    m.add_class::<variant::VariantFileFetchIter>()?;
    m.add_class::<variant::VariantRecord>()?;
    m.add_class::<variant::VariantSamples>()?;
    m.add_class::<variant::VariantSamplesIter>()?;
    m.add_class::<variant::VariantSample>()?;
    m.add_class::<variant::VariantHeader>()?;
    m.add_class::<variant::VariantContigs>()?;
    m.add_class::<variant::VariantContigsIter>()?;
    m.add_class::<variant::VariantContig>()?;
    m.add_class::<variant::VariantInfoDefs>()?;
    m.add_class::<variant::VariantInfoDefsIter>()?;
    m.add_class::<variant::VariantFormatDefs>()?;
    m.add_class::<variant::VariantFormatDefsIter>()?;
    m.add_class::<variant::VariantFieldDef>()?;
    Ok(())
}
