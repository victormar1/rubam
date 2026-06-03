//! Pure-Rust ports of samtools subcommands.
//!
//! Each public function (sort, index, view, merge, flagstat, idxstats,
//! calmd, faidx) is also re-exposed via the `rubam.tools` Python namespace
//! and the `rubam samtools` shadow CLI binary. Implementations land in
//! Phase B tasks B2..B8.

pub mod bcftools;
pub mod calmd;
pub mod faidx;
pub mod flagstat;
pub mod idxstats;
pub mod index;
pub mod merge;
pub mod sort;
pub mod view;
