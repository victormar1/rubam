//! `rubam-synth-bam` — pure-Rust synthetic BAM generator for benchmarking.
//!
//! Produces a coordinate-sorted, indexed BAM with single-end reads uniformly
//! distributed across one chromosome. No external dependency: Cargo + rustc
//! is enough to obtain a reproducible BAM of any size on Windows / Linux /
//! macOS.
//!
//! Usage:
//!     rubam-synth-bam --output sample.bam \
//!         --chrom chr20 --length 64444167 --coverage 30 \
//!         --read-length 150 --seed 42

use std::env;
use std::fs::File;
use std::io::BufWriter;
use std::num::NonZeroUsize;

use noodles::bam;
use noodles::core::Position;
use noodles::sam::{
    self,
    alignment::io::Write as _,
    alignment::record::cigar::op::Kind,
    alignment::record::cigar::Op,
    alignment::record::Flags,
    alignment::record::MappingQuality,
    alignment::record_buf::{Cigar, QualityScores, Sequence},
    alignment::RecordBuf,
    header::record::value::map::header::sort_order::COORDINATE,
    header::record::value::map::header::tag::SORT_ORDER,
    header::record::value::map::{self as map_kind, ReferenceSequence},
    header::record::value::Map,
};

/// xorshift64* — small, fast, good-enough PRNG for synthetic data.
struct Xor(u64);
impl Xor {
    fn new(seed: u64) -> Self {
        let s = seed.wrapping_mul(0x9E3779B97F4A7C15) | 1;
        Self(s)
    }
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }
    fn range(&mut self, n: u64) -> u64 {
        if n == 0 {
            0
        } else {
            self.next_u64() % n
        }
    }
}

#[derive(Debug)]
struct Args {
    output: String,
    chrom: String,
    chrom_len: u64,
    coverage: u32,
    read_len: usize,
    seed: u64,
    /// Probability that a read is spliced (CIGAR contains a single `N` op).
    /// 0.0 = no splicing (default, plain WGS-like). 0.05 = ~5 % spliced reads
    /// — a realistic figure for a deeply intronic RNA-seq region.
    spliced_frac: f64,
    /// Min/max intron length when `spliced_frac > 0`, in bp.
    intron_min: u64,
    intron_max: u64,
}

impl Args {
    fn parse() -> Result<Self, String> {
        let mut output: Option<String> = None;
        let mut chrom = "chr20".to_string();
        let mut chrom_len: u64 = 1_000_000;
        let mut coverage: u32 = 30;
        let mut read_len: usize = 150;
        let mut seed: u64 = 42;
        let mut spliced_frac: f64 = 0.0;
        let mut intron_min: u64 = 1_000;
        let mut intron_max: u64 = 10_000;

        let argv: Vec<String> = env::args().skip(1).collect();
        let mut i = 0;
        while i < argv.len() {
            let next = || {
                argv.get(i + 1)
                    .cloned()
                    .ok_or_else(|| format!("missing value for {}", argv[i]))
            };
            match argv[i].as_str() {
                "--output" | "-o" => {
                    output = Some(next()?);
                    i += 2;
                }
                "--chrom" => {
                    chrom = next()?;
                    i += 2;
                }
                "--length" => {
                    chrom_len = next()?.parse().map_err(|e| format!("--length: {e}"))?;
                    i += 2;
                }
                "--coverage" => {
                    coverage = next()?.parse().map_err(|e| format!("--coverage: {e}"))?;
                    i += 2;
                }
                "--read-length" => {
                    read_len = next()?.parse().map_err(|e| format!("--read-length: {e}"))?;
                    i += 2;
                }
                "--seed" => {
                    seed = next()?.parse().map_err(|e| format!("--seed: {e}"))?;
                    i += 2;
                }
                "--spliced" => {
                    spliced_frac = next()?.parse().map_err(|e| format!("--spliced: {e}"))?;
                    i += 2;
                }
                "--intron-min" => {
                    intron_min = next()?.parse().map_err(|e| format!("--intron-min: {e}"))?;
                    i += 2;
                }
                "--intron-max" => {
                    intron_max = next()?.parse().map_err(|e| format!("--intron-max: {e}"))?;
                    i += 2;
                }
                "--help" | "-h" => {
                    eprintln!(
                        "rubam-synth-bam — synthetic BAM generator\n\n\
                         Usage: rubam-synth-bam --output PATH --chrom NAME --length BP \\\n\
                         \t[--coverage 30] [--read-length 150] [--seed 42] \\\n\
                         \t[--spliced 0.05] [--intron-min 1000] [--intron-max 10000]\n\n\
                         When --spliced FRAC > 0, that fraction of reads are emitted\n\
                         with a CIGAR of the form `aM bN cM` (single intron skip),\n\
                         simulating spliced RNA-seq alignments. The cumulative `aM + cM`\n\
                         equals --read-length so the query stays consistent."
                    );
                    std::process::exit(0);
                }
                other => return Err(format!("unknown arg: {other}")),
            }
        }
        if !(0.0..=1.0).contains(&spliced_frac) {
            return Err("--spliced must be in [0, 1]".into());
        }
        if intron_max < intron_min {
            return Err("--intron-max must be >= --intron-min".into());
        }
        Ok(Self {
            output: output.ok_or_else(|| "--output is required".to_string())?,
            chrom,
            chrom_len,
            coverage,
            read_len,
            seed,
            spliced_frac,
            intron_min,
            intron_max,
        })
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse().map_err(|e| {
        eprintln!("error: {e}");
        e
    })?;

    let n_reads = (args.chrom_len * args.coverage as u64) / args.read_len as u64;
    eprintln!(
        "[synth] generating {} single-end reads ({}x cov over {} bp, read length {})",
        n_reads, args.coverage, args.chrom_len, args.read_len,
    );

    // Build header.
    let chrom_len_nz = NonZeroUsize::new(args.chrom_len as usize).ok_or("--length must be > 0")?;
    let header = sam::Header::builder()
        .set_header(
            Map::<map_kind::Header>::builder()
                .insert(SORT_ORDER, COORDINATE)
                .build()?,
        )
        .add_reference_sequence(
            args.chrom.as_str(),
            Map::<ReferenceSequence>::new(chrom_len_nz),
        )
        .build();

    // Generate sorted start positions. For spliced reads we also need a
    // bigger reference span (read_len + intron_max), so reserve room.
    let span_per_read = if args.spliced_frac > 0.0 {
        args.read_len as u64 + args.intron_max
    } else {
        args.read_len as u64
    };
    let max_start = args.chrom_len.saturating_sub(span_per_read) + 1;
    if max_start < 1 {
        return Err(format!(
            "read_length+intron ({}) > chrom length ({})",
            span_per_read, args.chrom_len
        )
        .into());
    }
    let mut rng = Xor::new(args.seed);
    let mut starts: Vec<u64> = (0..n_reads).map(|_| rng.range(max_start) + 1).collect();
    starts.sort_unstable();

    // Open BAM writer.
    let bam_path = args.output.clone();
    let writer_file = BufWriter::new(File::create(&bam_path)?);
    let mut bam_writer = bam::io::Writer::new(writer_file);
    bam_writer.write_header(&header)?;

    let plain_cigar: Cigar = [Op::new(Kind::Match, args.read_len)].into_iter().collect();
    let mapq = MappingQuality::new(60).expect("60 is a valid MAPQ");
    let qual_template = QualityScores::from(vec![30u8; args.read_len]);

    let mut record = RecordBuf::default();
    let mut n_spliced: u64 = 0;
    let intron_range = args.intron_max - args.intron_min + 1;

    for (idx, start) in starts.iter().enumerate() {
        let mut seq_bytes = vec![0u8; args.read_len];
        for j in 0..args.read_len {
            seq_bytes[j] = match rng.range(4) {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };
        }
        // Decide whether this read is spliced. Use the PRNG for determinism.
        // We compare a uniform u64 against threshold = spliced_frac * 2^64.
        let r01 = (rng.next_u64() as f64) / (u64::MAX as f64);
        let cigar = if r01 < args.spliced_frac {
            // Pick a split point (>= 5 bp on each side) and an intron length.
            let min_arm = 5usize.min(args.read_len.saturating_sub(5));
            let max_arm = args.read_len.saturating_sub(5);
            let split = if max_arm > min_arm {
                min_arm + (rng.range((max_arm - min_arm) as u64) as usize)
            } else {
                args.read_len / 2
            };
            let intron_len = args.intron_min + rng.range(intron_range);
            n_spliced += 1;
            [
                Op::new(Kind::Match, split),
                Op::new(Kind::Skip, intron_len as usize),
                Op::new(Kind::Match, args.read_len - split),
            ]
            .into_iter()
            .collect::<Cigar>()
        } else {
            plain_cigar.clone()
        };

        *record.name_mut() = Some(format!("synth-{idx}").into());
        *record.flags_mut() = Flags::default();
        *record.reference_sequence_id_mut() = Some(0);
        *record.alignment_start_mut() = Position::new(*start as usize);
        *record.mapping_quality_mut() = Some(mapq);
        *record.cigar_mut() = cigar;
        *record.sequence_mut() = Sequence::from(seq_bytes);
        *record.quality_scores_mut() = qual_template.clone();
        *record.template_length_mut() = 0;
        *record.mate_reference_sequence_id_mut() = None;
        *record.mate_alignment_start_mut() = None;

        bam_writer.write_alignment_record(&header, &record)?;

        if idx > 0 && idx.is_multiple_of(1_000_000) {
            eprintln!("[synth]   {} / {}", idx, n_reads);
        }
    }
    if args.spliced_frac > 0.0 {
        eprintln!(
            "[synth] spliced reads: {} / {} ({:.2}%)",
            n_spliced,
            n_reads,
            100.0 * n_spliced as f64 / n_reads as f64
        );
    }
    bam_writer.finish(&header)?;
    drop(bam_writer);
    eprintln!("[synth] wrote {n_reads} records to {bam_path}");

    // Build index.
    eprintln!("[synth] indexing -> {bam_path}.bai");
    let index = bam::fs::index(&bam_path)?;
    let bai_path = format!("{bam_path}.bai");
    let bai_file = File::create(&bai_path)?;
    let mut bai_writer = bam::bai::io::Writer::new(bai_file);
    bai_writer.write_index(&index)?;
    eprintln!("[synth] done.");
    Ok(())
}
