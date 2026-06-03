//! Aux — typed enum mirroring rust_htslib::bam::record::Aux.
//!
//! Bridge from noodles' Value (untyped Cow-based) to a typed enum keeps
//! HARMOS's pattern-matching code drop-in.

use thiserror::Error;

/// A typed BAM auxiliary tag value.
///
/// Variant names and shapes match `rust_htslib::bam::record::Aux` so
/// HARMOS's `record.aux(b"SA") -> Aux::String(s)` pattern works
/// unchanged. Borrowed variants tie their lifetime to the parsed
/// record's backing buffer.
#[derive(Clone, Debug, PartialEq)]
pub enum Aux<'a> {
    Char(u8),
    I8(i8),
    U8(u8),
    I16(i16),
    U16(u16),
    I32(i32),
    U32(u32),
    Float(f32),
    Double(f64),
    String(&'a str),
    HexByteArray(&'a str),
    ArrayI8(&'a [i8]),
    ArrayU8(&'a [u8]),
    ArrayI16(&'a [i16]),
    ArrayU16(&'a [u16]),
    ArrayI32(&'a [i32]),
    ArrayU32(&'a [u32]),
    ArrayFloat(&'a [f32]),
}

/// Errors produced while looking up or decoding an aux tag.
#[derive(Debug, Error)]
pub enum AuxError {
    #[error("tag name must be 2 ASCII bytes, got {0} bytes")]
    BadTagLength(usize),
    #[error("tag {0:?} not present in record")]
    NotFound(String),
    #[error("tag {0:?} value parse error: {1}")]
    Parse(String, String),
    #[error("unsupported aux variant for tag {0:?}: {1}")]
    Unsupported(String, String),
}

/// Convert a `noodles` field Value into our typed `Aux`.
///
/// Scalar/string variants borrow directly from the noodles value, which itself
/// borrows from the parent record. No allocation, no arena, no leak.
///
/// Array variants will require an owned-vector arena (see `AuxArena` below)
/// when they ship in v0.4; until then they return `AuxError::Unsupported`.
///
/// **v0.3.2 (Wave 3 of the major revision)**: the previous signature took
/// `_arena: &'a AuxArena` and forced callers to leak a `Box<AuxArena>` to
/// obtain a `'static` lifetime. That added ~24 bytes of heap garbage per
/// `aux()` call and was the reviewer's M1 finding. The arena parameter is
/// gone now; lifetimes flow purely from the noodles value.
pub(crate) fn aux_from_noodles<'a>(
    tag_name: [u8; 2],
    value: noodles::sam::alignment::record::data::field::Value<'a>,
) -> Result<Aux<'a>, AuxError> {
    use noodles::sam::alignment::record::data::field::Value as V;
    let tag_str = std::str::from_utf8(&tag_name).unwrap_or("??").to_string();
    Ok(match value {
        V::Character(c) => Aux::Char(c),
        V::Int8(n) => Aux::I8(n),
        V::UInt8(n) => Aux::U8(n),
        V::Int16(n) => Aux::I16(n),
        V::UInt16(n) => Aux::U16(n),
        V::Int32(n) => Aux::I32(n),
        V::UInt32(n) => Aux::U32(n),
        V::Float(n) => Aux::Float(n),
        V::String(s) => Aux::String(
            std::str::from_utf8(s.as_ref())
                .map_err(|e| AuxError::Parse(tag_str.clone(), e.to_string()))?,
        ),
        V::Hex(s) => Aux::HexByteArray(
            std::str::from_utf8(s.as_ref())
                .map_err(|e| AuxError::Parse(tag_str.clone(), e.to_string()))?,
        ),
        V::Array(_arr) => {
            return Err(AuxError::Unsupported(
                tag_str,
                "array tags land in v0.4 with the proper owned-vector arena. \
                 HARMOS does not currently use array tags in the SA-tag path."
                    .to_string(),
            ));
        }
    })
}

/// Reserved storage for owned-vector array tags. Currently zero-sized; the
/// type is preserved across the v0.3.2 refactor for forward-compatibility so
/// that v0.4 can introduce array-tag support without breaking the public API
/// surface of `aux_data`.
///
/// **No allocation, no leak**: an `AuxArena::default()` is a zero-byte struct.
#[derive(Default, Debug)]
pub struct AuxArena {
    _placeholder: (),
}
