//! Pure-Rust ports of bcftools subcommands.
//!
//! Mirrors `src/tools/` (samtools side). Each subcommand file exposes a
//! `_native` Rust function (`io::Result`) and a pyfunction wrapper. The
//! `rubam-bcftools` shadow CLI binary in `src/bin/bcftools.rs` calls the
//! native fns. Phase B of `paper/PLAN_v0.3.md` implements them.

pub mod concat;
pub mod index;
pub mod norm;
pub mod query;
pub mod sort;
pub mod stats;
pub mod view;
