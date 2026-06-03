//! AlignmentFile — opens a BAM and iterates records.
//!
//! Linear iteration only in v0.2.1 (HARMOS doesn't use IndexedReader).
//! The file is opened in streaming mode regardless of whether `.bai` /
//! `.csi` are present, mirroring `rust_htslib::bam::Reader::from_path`.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use noodles::bam;

use super::aligned_segment::AlignedSegment;
use super::error::{Error, Result};
use super::header::Header;

/// Streaming BAM reader. `records()` yields `AlignedSegment`s in file order.
pub struct AlignmentFile {
    reader: bam::io::Reader<noodles::bgzf::io::Reader<Box<dyn BufRead + Send>>>,
    header: Header,
}

impl AlignmentFile {
    /// Open a BAM by path. Auto-detects BGZF compression. Does **not** require
    /// a `.bai` / `.csi` index — HARMOS iterates linearly, no random access.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let f: Box<dyn BufRead + Send> = Box::new(BufReader::new(File::open(path)?));
        // bam::io::Reader::new wraps the inner reader in bgzf::io::Reader internally.
        let mut reader = bam::io::Reader::new(f);
        let h = crate::common::read_bam_header_tolerant(reader.get_mut()).map_err(Error::Io)?;
        Ok(Self {
            reader,
            header: Header::from_noodles(h),
        })
    }

    /// The parsed header. Exposes `tid2name`, `target_count`, etc.
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Iterator over alignment records. Each item is `Result<AlignedSegment>`
    /// so HARMOS can choose to ignore parse errors per-record (`.flatten()`)
    /// or fail fast (`.collect::<Result<Vec<_>, _>>()`).
    pub fn records(&mut self) -> Records<'_> {
        Records {
            reader: &mut self.reader,
            header: &self.header,
        }
    }
}

/// Iterator returned by [`AlignmentFile::records`].
pub struct Records<'a> {
    reader: &'a mut bam::io::Reader<noodles::bgzf::io::Reader<Box<dyn BufRead + Send>>>,
    header: &'a Header,
}

impl<'a> Iterator for Records<'a> {
    type Item = Result<AlignedSegment>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = bam::Record::default();
        match self.reader.read_record(&mut record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(AlignedSegment::new(record, self.header.clone()))),
            Err(e) => Some(Err(Error::Io(e))),
        }
    }
}
