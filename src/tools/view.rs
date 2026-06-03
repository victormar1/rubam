//! samtools view — region/flag/MAPQ filter and optional BAM output.

use std::fs::File;
use std::io::BufWriter;

use noodles::bam;
use noodles::bgzf;
use noodles::core::Region;
#[cfg(feature = "python")]
use pyo3::exceptions::{PyIOError, PyValueError};
#[cfg(feature = "python")]
use pyo3::prelude::*;

use crate::common::{open_indexed, open_streaming, read_header_indexed, read_header_streaming};

/// samtools view equivalent.
///
/// If `count_only=True`, returns the matching record count and ignores `output`.
/// Otherwise, if `output` is set, writes a filtered BAM and still returns the count.
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (input, *, region = None, output = None,
                    min_mapq = 0, flag_required = 0, flag_filtered = 0,
                    count_only = false))]
pub fn view(
    input: &str,
    region: Option<&str>,
    output: Option<&str>,
    min_mapq: u8,
    flag_required: u16,
    flag_filtered: u16,
    count_only: bool,
) -> PyResult<u64> {
    view_native(
        input,
        region,
        output,
        min_mapq,
        flag_required,
        flag_filtered,
        count_only,
    )
    .map_err(|e| match e.kind() {
        std::io::ErrorKind::InvalidInput => PyValueError::new_err(e.to_string()),
        _ => PyIOError::new_err(format!("view: {e}")),
    })
}

/// Pure-Rust entry point used by the rubam-samtools shadow CLI in Phase C.
pub fn view_native(
    input: &str,
    region: Option<&str>,
    output: Option<&str>,
    min_mapq: u8,
    flag_required: u16,
    flag_filtered: u16,
    count_only: bool,
) -> std::io::Result<u64> {
    fn pass(flags: u16, mq: u8, fr: u16, ff: u16, mq_min: u8) -> bool {
        if (flags & fr) != fr {
            return false;
        }
        if (flags & ff) != 0 {
            return false;
        }
        if mq < mq_min {
            return false;
        }
        true
    }

    type BamWriter = bam::io::Writer<bgzf::io::Writer<BufWriter<File>>>;

    fn make_writer(out: &str, header: &noodles::sam::Header) -> std::io::Result<BamWriter> {
        let f = File::create(out)?;
        let mut w = bam::io::Writer::new(BufWriter::new(f));
        w.write_header(header)?;
        Ok(w)
    }

    if let Some(reg_str) = region {
        // Indexed path.
        let mut reader = open_indexed(input)?;
        let header = read_header_indexed(&mut reader)?;
        let parsed: Region = reg_str
            .parse()
            .map_err(|e: noodles::core::region::ParseError| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    format!("parse region {reg_str:?}: {e}"),
                )
            })?;
        let mut writer_opt: Option<BamWriter> = if !count_only {
            output.map(|out| make_writer(out, &header)).transpose()?
        } else {
            None
        };
        let mut count: u64 = 0;
        let query = reader.query(&header, &parsed)?;
        for r in query.records() {
            let r = r?;
            let flags = r.flags().bits();
            let mq = r.mapping_quality().map(|q| q.get()).unwrap_or(0);
            if !pass(flags, mq, flag_required, flag_filtered, min_mapq) {
                continue;
            }
            count += 1;
            if let Some(w) = writer_opt.as_mut() {
                w.write_record(&header, &r)?;
            }
        }
        if let Some(mut w) = writer_opt {
            w.try_finish()?;
        }
        return Ok(count);
    }

    // Streaming path.
    let mut reader = open_streaming(input)?;
    let header = read_header_streaming(&mut reader)?;
    let mut writer_opt: Option<BamWriter> = if !count_only {
        output.map(|out| make_writer(out, &header)).transpose()?
    } else {
        None
    };
    let mut count: u64 = 0;
    let mut buf = bam::Record::default();
    loop {
        match reader.read_record(&mut buf)? {
            0 => break,
            _ => {
                let flags = buf.flags().bits();
                let mq = buf.mapping_quality().map(|q| q.get()).unwrap_or(0);
                if !pass(flags, mq, flag_required, flag_filtered, min_mapq) {
                    continue;
                }
                count += 1;
                if let Some(w) = writer_opt.as_mut() {
                    w.write_record(&header, &buf)?;
                }
            }
        }
    }
    if let Some(mut w) = writer_opt {
        w.try_finish()?;
    }
    Ok(count)
}
