//! Shared helpers: BAM opening, header parsing, default flag mask.
//!
//! Pure-Rust: returns `std::io::Result<_>`. The pyo3 wrappers translate
//! to `PyIOError` / `PyValueError` at the FFI boundary.

use std::fs::File;
use std::io::{self, Read};
use std::num::NonZero;

use noodles::bam;
use noodles::bgzf;
use noodles::sam::header::{
    record::value::{map::ReferenceSequence, Map},
    Parser, ReferenceSequences,
};

/// Default SAM-flag mask: `UNMAP (0x4) | SECONDARY (0x100) | QCFAIL (0x200) | DUP (0x400)`.
pub const FLAG_FILTER_DEFAULT: u16 = 0x4 | 0x100 | 0x200 | 0x400;

/// BAM magic number (`BAM\1`).
const BAM_MAGIC: &[u8; 4] = b"BAM\x01";

pub type IndexedBamReader = bam::io::IndexedReader<bgzf::io::Reader<File>>;
pub type StreamingBamReader = bam::io::Reader<bgzf::io::Reader<File>>;

pub fn open_indexed(path: &str) -> io::Result<IndexedBamReader> {
    bam::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .map_err(|e| {
            io::Error::new(
                io::ErrorKind::Other,
                format!("failed to open indexed BAM at {path}: {e}"),
            )
        })
}

pub fn open_streaming(path: &str) -> io::Result<StreamingBamReader> {
    bam::io::reader::Builder.build_from_path(path).map_err(|e| {
        io::Error::new(
            io::ErrorKind::Other,
            format!("failed to open BAM at {path}: {e}"),
        )
    })
}

pub fn read_header_indexed(reader: &mut IndexedBamReader) -> io::Result<noodles::sam::Header> {
    read_bam_header_tolerant(reader.get_mut())
        .map_err(|e| io::Error::new(io::ErrorKind::Other, format!("failed to read header: {e}")))
}

pub fn read_header_streaming(reader: &mut StreamingBamReader) -> io::Result<noodles::sam::Header> {
    read_bam_header_tolerant(reader.get_mut())
        .map_err(|e| io::Error::new(io::ErrorKind::Other, format!("failed to read header: {e}")))
}

fn read_i32_le<R: Read>(r: &mut R) -> io::Result<i32> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf)?;
    Ok(i32::from_le_bytes(buf))
}

/// Coerce an arbitrary `@HD` `VN` value into the strict `MAJOR.MINOR` form
/// that noodles requires. `1.6.0` → `1.6`, `1` → `1.0`, garbage → `1.6`.
fn coerce_version(v: &str) -> String {
    let parts: Vec<&str> = v.split('.').collect();
    if parts.len() >= 2 && parts[0].parse::<u32>().is_ok() && parts[1].parse::<u32>().is_ok() {
        format!("{}.{}", parts[0], parts[1])
    } else if parts.len() == 1 && !parts[0].is_empty() && parts[0].parse::<u32>().is_ok() {
        format!("{}.0", parts[0])
    } else {
        "1.6".to_string()
    }
}

/// Rewrite an `@HD` line so it always carries a strict `VN:MAJOR.MINOR` field
/// (htslib/pysam tolerate a missing or malformed version; noodles does not).
fn sanitize_hd_line(line: &[u8]) -> Vec<u8> {
    let s = String::from_utf8_lossy(line);
    let mut fields: Vec<String> = s.split('\t').map(|f| f.to_string()).collect();
    match fields.iter().position(|f| f.starts_with("VN:")) {
        Some(i) => {
            let coerced = coerce_version(&fields[i][3..]);
            fields[i] = format!("VN:{coerced}");
        }
        None => {
            // Insert a default version right after the `@HD` tag.
            let at = fields.len().min(1);
            fields.insert(at, "VN:1.6".to_string());
        }
    }
    fields.join("\t").into_bytes()
}

/// Iterate the SAM-text header lines: split on `\n`, strip a trailing `\r`,
/// and drop anything that is not a real header record (`@`-prefixed). This
/// also transparently skips any NUL padding htslib appends to the text block.
fn header_lines(text: &[u8]) -> impl Iterator<Item = &[u8]> {
    text.split(|&b| b == b'\n').filter_map(|line| {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.first() == Some(&b'@') {
            Some(line)
        } else {
            None
        }
    })
}

/// Tolerant BAM header reader.
///
/// Reads the raw BAM header structure (magic, SAM text, binary reference
/// dictionary) directly from `reader`, leaving the stream positioned at the
/// first alignment record — exactly like `noodles`' own `read_header`, but
/// without its strict SAM-header parsing.
///
/// htslib/pysam accept real-world headers that noodles' strict parser rejects
/// with `invalid record`: an `@HD` line with no `VN` field or a multi-part
/// version (`VN:1.6.0`), duplicate `@PG`/`@RG`/`@SQ` IDs (common in re-run
/// GATK/Picard pipelines), etc. To replace pysam on full hg38 BAMs we must be
/// equally tolerant.
///
/// Strategy: try the strict parser first (zero behaviour change for
/// well-formed files, preserving every `@SQ`/`@RG`/`@PG` tag); on failure,
/// fall back to a best-effort parse that sanitizes `@HD`, skips unparseable
/// lines, and rebuilds the reference dictionary from the authoritative binary
/// list that BAM always carries.
pub fn read_bam_header_tolerant<R: Read>(reader: &mut R) -> io::Result<noodles::sam::Header> {
    // ---- magic ----
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;
    if &magic != BAM_MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BAM magic number",
        ));
    }

    // ---- SAM text header ----
    let l_text = read_i32_le(reader)?;
    let mut text = vec![0u8; l_text.max(0) as usize];
    reader.read_exact(&mut text)?;

    // ---- binary reference dictionary (authoritative) ----
    let n_ref = read_i32_le(reader)?;
    let mut bin_refs = ReferenceSequences::default();
    for _ in 0..n_ref.max(0) {
        let l_name = read_i32_le(reader)?.max(0) as usize;
        let mut name = vec![0u8; l_name];
        reader.read_exact(&mut name)?;
        while name.last() == Some(&0) {
            name.pop();
        }
        let l_ref = read_i32_le(reader)?;
        let len = NonZero::new(l_ref.max(1) as usize).unwrap();
        bin_refs.insert(name.into(), Map::<ReferenceSequence>::new(len));
    }

    // ---- fast path: strict parse (no regression for valid files) ----
    let mut strict = Parser::default();
    let mut strict_ok = true;
    for line in header_lines(&text) {
        if strict.parse_partial(line).is_err() {
            strict_ok = false;
            break;
        }
    }
    if strict_ok {
        let mut header = strict.finish();
        if header.reference_sequences().is_empty() {
            *header.reference_sequences_mut() = bin_refs;
        }
        return Ok(header);
    }

    // ---- tolerant fallback ----
    let mut parser = Parser::default();
    for line in header_lines(&text) {
        if parser.parse_partial(line).is_err() && line.starts_with(b"@HD") {
            // The only structurally recoverable failure is a bad/missing @HD
            // version; everything else (dup IDs, exotic lines) is simply
            // dropped — htslib keeps such records but we only need a usable,
            // non-crashing header.
            let _ = parser.parse_partial(&sanitize_hd_line(line));
        }
    }
    let mut header = parser.finish();
    // The binary dictionary is the source of truth for contig names/lengths.
    *header.reference_sequences_mut() = bin_refs;
    Ok(header)
}

pub fn validate_chrom(header: &noodles::sam::Header, chrom: &str) -> io::Result<()> {
    if header.reference_sequences().contains_key(chrom.as_bytes()) {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("chromosome {chrom} not found in BAM header"),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Assemble the raw (uncompressed) BAM header structure for the given SAM
    /// text and `(name, length)` reference list.
    fn raw_bam_header(sam_text: &str, refs: &[(&str, i32)]) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(BAM_MAGIC);
        out.extend_from_slice(&(sam_text.len() as i32).to_le_bytes());
        out.extend_from_slice(sam_text.as_bytes());
        out.extend_from_slice(&(refs.len() as i32).to_le_bytes());
        for (name, len) in refs {
            let mut nm = name.as_bytes().to_vec();
            nm.push(0);
            out.extend_from_slice(&(nm.len() as i32).to_le_bytes());
            out.extend_from_slice(&nm);
            out.extend_from_slice(&len.to_le_bytes());
        }
        out
    }

    const REFS: &[(&str, i32)] = &[
        ("chr1", 248_956_422),
        ("chr2", 242_193_529),
        ("chrM", 16_569),
    ];

    fn sq_block() -> String {
        REFS.iter()
            .map(|(n, l)| format!("@SQ\tSN:{n}\tLN:{l}\n"))
            .collect()
    }

    fn assert_refs_ok(sam_text: &str) {
        let bytes = raw_bam_header(sam_text, REFS);
        let header = read_bam_header_tolerant(&mut &bytes[..]).expect("tolerant header parse");
        let names: Vec<String> = header
            .reference_sequences()
            .keys()
            .map(|k| String::from_utf8_lossy(k).into_owned())
            .collect();
        assert_eq!(names, vec!["chr1", "chr2", "chrM"]);
    }

    #[test]
    fn tolerant_hd_without_version() {
        assert_refs_ok(&format!("@HD\tSO:coordinate\n{}", sq_block()));
    }

    #[test]
    fn tolerant_multipart_version() {
        assert_refs_ok(&format!("@HD\tVN:1.6.0\tSO:coordinate\n{}", sq_block()));
    }

    #[test]
    fn tolerant_duplicate_program_id() {
        let txt = format!(
            "@HD\tVN:1.6\n{}@PG\tID:bwa\tCL:a\n@PG\tID:bwa\tCL:b\n",
            sq_block()
        );
        assert_refs_ok(&txt);
    }

    #[test]
    fn well_formed_header_preserves_metadata() {
        let txt = format!(
            "@HD\tVN:1.6\tSO:coordinate\n{}@RG\tID:s1\tSM:NA12878\n@PG\tID:bwa\tPN:bwa\n",
            sq_block()
        );
        let bytes = raw_bam_header(&txt, REFS);
        let header = read_bam_header_tolerant(&mut &bytes[..]).unwrap();
        // The strict fast path must run for valid headers, keeping @RG/@PG.
        assert_eq!(header.read_groups().len(), 1);
        assert_eq!(header.programs().as_ref().len(), 1);
        assert!(header.header().is_some());
    }

    #[test]
    fn coerce_version_cases() {
        assert_eq!(coerce_version("1.6"), "1.6");
        assert_eq!(coerce_version("1.6.0"), "1.6");
        assert_eq!(coerce_version("1"), "1.0");
        assert_eq!(coerce_version("garbage"), "1.6");
        assert_eq!(coerce_version(""), "1.6");
    }
}
