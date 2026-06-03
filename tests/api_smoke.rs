//! Smoke test: the public api::* re-exports compile.
//!
//! Each type is filled in by tasks A2..A6. This test only checks that the
//! imports resolve and that downstream Rust crates can `use rubam::api::*`.

use rubam::api::{AlignedSegment, AlignmentFile, Aux, Cigar, Header};

#[test]
fn types_compile() {
    let _ = std::any::TypeId::of::<AlignmentFile>();
    let _ = std::any::TypeId::of::<AlignedSegment>();
    let _ = std::any::TypeId::of::<Header>();
    let _ = std::any::TypeId::of::<Cigar>();
    // Aux carries a lifetime, so use a concrete instantiation for the type id.
    let _ = std::any::TypeId::of::<Aux<'static>>();
}
