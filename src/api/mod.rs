//! Pure-Rust public API. The pyo3 wrappers in `src/alignment.rs` etc. are
//! built on top of these. External Rust crates (HARMOS, etc.) depend on
//! these types directly without pulling in pyo3.

pub mod aligned_segment;
pub mod alignment_file;
pub mod aux_data;
pub mod cigar;
pub mod error;
pub mod header;

pub use aligned_segment::AlignedSegment;
pub use alignment_file::AlignmentFile;
pub use aux_data::{Aux, AuxError};
pub use cigar::Cigar;
pub use error::{Error, Result};
pub use header::Header;
