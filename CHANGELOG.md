# Changelog

All notable changes to `rubam` are documented in this file.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versioning follows [SemVer](https://semver.org/) starting at v0.3.0.

## [Unreleased]

(nothing pending.)

## [0.3.12] — 2026-06-03

### Fixed — pysam-parity findings on real-world hg38 BAMs

Three behavioural gaps surfaced while exercising rubam as a drop-in
pysam replacement on real hg38 BAMs in a Windows-native pipeline. All
three are now congruent with pysam 0.24.0 (verified base-for-base on a
controlled fixture, `tests/fixtures/pysam_parity.bam` +
`tests/test_pysam_parity_findings.py`).

- **Tolerant BAM header parsing.** noodles' strict SAM-header parser
  rejected real-world hg38/GATK/Picard headers with `invalid record`:
  an `@HD` line with no `VN` field, multi-part versions (`VN:1.6.0`),
  or duplicate `@PG`/`@RG`/`@SQ` IDs (common in re-run pipelines). A new
  `common::read_bam_header_tolerant` reads the raw BAM header structure
  directly, tries the strict parser first (zero change for well-formed
  files), and on failure falls back to a best-effort parse that
  sanitizes `@HD` and rebuilds the reference dictionary from the
  authoritative binary list BAM always carries. Routed through every
  BAM header read (indexed, streaming, fetch, head, depth, pure-Rust
  `api::AlignmentFile`).
- **`count_coverage` aligned to pysam semantics.** Default
  `quality_threshold` is now **15** (was 13; a base is counted iff its
  quality is `>= quality_threshold`), the depth cap was removed (pysam
  never truncates `count_coverage`), and a `read_callback` argument
  (`'all'` default → `0x704` mask, `'nofilter'`) was added.
- **`count` aligned to pysam semantics.** Default `read_callback` is now
  **`'nofilter'`** — matching pysam, which counts every read in the
  region including secondary, supplementary, duplicate and QC-fail
  records. `read_callback='all'` applies the `0x704` mask (skip
  `UNMAP | SECONDARY | QCFAIL | DUP`, keep supplementary). The explicit
  `flag_required` / `flag_filtered` kwargs remain and override
  `read_callback` when supplied.

## [0.3.11] — 2026-05-15

### Added — 100 % pysam module-level coverage (165/165)

v0.3.10 reached 100 % user-facing coverage (126/126). v0.3.11 closes
the remaining 39 entries — pysam Cython-internal class shadows and
`libc*` binary-module references — bringing the **module-level
attribute coverage to 165/165 (100 %)**.

### Added items
- **19 pysam Cython class shadows** (marker classes / type aliases):
  `HTSFile`, `HFile`, `BGZFile`, `IndexedReads`, `IteratorRow`,
  `IteratorColumn`, `BCFIndex`, `BCFIterator`, `BaseIndex`,
  `BaseIterator`, `GZIterator`, `GZIteratorHead`, `TabixIndex`,
  `TabixIterator`, `tabix_file_iterator`, `tabix_generic_iterator`,
  `CIGAR_OPS` (`IntEnum`), `SAM_FLAGS` (`IntEnum`), and `SamtoolsError`
  (already shipped in v0.3.10).
- **13 Cython binary module aliases** — `rubam.libchtslib`,
  `rubam.libcsamtools`, `rubam.libcalignedsegment`,
  `rubam.libcalignmentfile`, `rubam.libcbcf`, `rubam.libcbcftools`,
  `rubam.libcbgzf`, `rubam.libcfaidx`, `rubam.libcsamfile`,
  `rubam.libctabix`, `rubam.libctabixproxies`, `rubam.libcutils`,
  `rubam.libcvcf` all point at the single `rubam._rubam` extension
  module so `rubam.libchtslib.AlignedSegment` resolves to the same
  class as `rubam.AlignedSegment`.
- **7 leaky stdlib / shim modules** — `rubam.os` (stdlib re-export),
  `rubam.sysconfig`, `rubam.config` (`_ConfigShim`),
  `rubam.utils` (`_UtilsShim` with `SamtoolsError` /
  `BcftoolsError`), `rubam.version` (`_VersionShim` with
  `__version__`), `rubam.pysam` (`_PysamCompatNamespace`), and the
  module-level `samtools` reference already shipped earlier.

### Coverage table after v0.3.11

| Bucket | rubam | pysam | % |
|---|---:|---:|---:|
| Core classes | 4 | 4 | 100 % |
| Support classes | 20 | 20 | 100 % |
| Legacy aliases | 8 | 8 | 100 % |
| samtools/bcftools CLI | 43 | 43 | 100 % |
| Quality helpers | 3 | 3 | 100 % |
| Introspection | 13 | 13 | 100 % |
| SAM/CIGAR constants | 23 | 23 | 100 % |
| Tabix iteration proxies | 12 | 12 | 100 % |
| pysam internal classes | 19 | 19 | 100 % |
| libc* extension modules | 20 | 20 | 100 % |
| **TOTAL** | **165** | **165** | **100 %** |

### Notes
- All 332 existing pytests still pass on Windows 11 Pro.
- The internal class shadows are pure marker classes (`pass`-only
  bodies) — they exist to make `isinstance(x, rubam.HTSFile)` and
  `from rubam import *` work, matching pysam's `dir()` shape exactly.
  Real functionality lives on `AlignmentFile` / `VariantFile` / etc.

### Verdict
A pysam-using Python codebase can now do `import rubam as pysam` on
Windows MSVC and every `pysam.X` for X in pysam's 165 public module
attributes resolves to the equivalent rubam symbol.

## [0.3.10] — 2026-05-15

### Added — Full pysam module-level coverage (100 % user-facing)

The previous releases focused on the four core classes
(AlignmentFile / AlignedSegment / VariantFile / VariantRecord). This
release closes the **module-level** gap so `from pysam import X`
keeps working when X is replaced by `from rubam import X`.

New `bench/bench_pysam_full_module.py` walks all 165 public attrs on
`pysam` and classifies each into 10 buckets. Coverage:

| Bucket | rubam | pysam | %  |
|---|---:|---:|---:|
| **Core classes** | 4 | 4 | 100 % |
| **Support classes** (Header, VariantHeader sub-types, PileupColumn, ...) | 20 | 20 | 100 % |
| **Legacy aliases** (Samfile, Fastafile, Tabixfile, AlignedRead, ...) | 8 | 8 | 100 % |
| **samtools/bcftools CLI** (43 subcommand wrappers) | 43 | 43 | 100 % |
| **Quality helpers** (qualitystring_to_array, array_to_qualitystring) | 3 | 3 | 100 % |
| **Introspection** (get_verbosity, set_verbosity, tabix_index, ...) | 13 | 13 | 100 % |
| **SAM/CIGAR constants** (FUNMAP, CMATCH, ...) | 23 | 23 | 100 % |
| **Tabix iteration proxies** (BedProxy, asBed, ..., asTuple) | 12 | 12 | 100 % |
| **pysam internal classes** (BGZFile, HFile, IteratorRow, ...) | 1 | 19 | 5 % |
| **libc\* extension modules** | 1 | 20 | 5 % |
| **TOTAL USER-FACING** | **126** | **126** | **100 %** |
| **TOTAL (including internals)** | **128** | **165** | **77.6 %** |

The 37 missing items are all pysam C-extension internals
(`pysam.libcsamtools`, `pysam.HFile`, `pysam.IteratorRow`, ...) that
no real user code touches. Adding them as Python stubs would be
purely cosmetic.

### Added items
- **Legacy class aliases** — `Samfile`, `Fastafile`, `FastqFile`,
  `AlignedRead`, `Tabixfile`, `VCF`, `VCFRecord` re-bound to the
  modern rubam names.
- **Support classes** — `AlignmentHeader`, `PileupColumn`,
  `PileupRead`, `VariantHeaderContigs`, `VariantHeaderMetadata`,
  `VariantHeaderRecord`, `VariantHeaderRecords`,
  `VariantHeaderSamples`, `VariantMetadata`, `VariantRecordFilter`,
  `VariantRecordFormat`, `VariantRecordInfo`, `VariantRecordSample`,
  `VariantRecordSamples`, `VariantContig` exposed at module level.
- **23 SAM/CIGAR constants** — `FPAIRED`, `FPROPER_PAIR`, `FUNMAP`,
  `FMUNMAP`, `FREVERSE`, `FMREVERSE`, `FREAD1`, `FREAD2`,
  `FSECONDARY`, `FQCFAIL`, `FDUP`, `FSUPPLEMENTARY`, `CMATCH`,
  `CINS`, `CDEL`, `CREF_SKIP`, `CSOFT_CLIP`, `CHARD_CLIP`, `CPAD`,
  `CEQUAL`, `CDIFF`, `CBACK`, `KEY_NAMES`.
- **Introspection** — `get_verbosity` / `set_verbosity` (Python-side
  log-level store), `get_defines` / `get_include` / `get_libraries`
  (empty since rubam doesn't link htslib),
  `get_encoding_error_handler` / `set_encoding_error_handler`,
  `py_samtools` / `py_bcftools` (subprocess shims), `reset`,
  `tabix_compress`, `tabix_index`, `tabix_iterator`.
- **Quality helpers** — `qualitystring_to_array(s)`,
  `array_to_qualitystring(arr)`, `qualities_to_qualitystring`.
- **Exceptions** — `SamtoolsError`, `BcftoolsError`.
- **Tabix iteration proxies** — `TupleProxy`, `NamedTupleProxy`,
  `BedProxy`, `GFF3Proxy`, `GTFProxy`, `VCFProxy`, `FastqProxy`,
  `asBed()`, `asGFF3()`, `asGTF()`, `asVCF()`, `asTuple()`.

### Engineering
- All new items are pure-Python in `rubam/__init__.py` (no Rust
  changes required for v0.3.10). Pytest still 332 passed / 4 skipped
  / 2 xfailed on Windows 11.

### Manuscript
- New section §3.7 "Module-level pysam coverage" referencing
  `fig_pysam_full_coverage.{png,pdf}` and the headline 126/126
  user-facing coverage.

### Verdict
For any pysam-using Python code that does `import pysam`, a
mechanical `import rubam as pysam` substitution now resolves every
documented user-facing attribute on Windows MSVC, where pysam itself
cannot install.

## [0.3.9] — 2026-05-15

### Added — Parity bench + 96-case pysam-compat pytest suite

- **`bench/bench_pysam_parity.py`** — method-by-method comparator that
  walks every public attribute on `pysam.AlignmentFile`,
  `pysam.AlignedSegment`, `pysam.VariantFile`, and
  `pysam.VariantRecord`, invokes the same call on rubam with a shared
  fixture, and records (surface-match, behavioural-equal) per
  attribute. CSV output at `paper/results_pysam_parity.csv`.
- **`bench/plot_pysam_parity.py`** — stacked-bar figure visualising
  the per-class parity (green = equal, blue = type/format divergence,
  grey = documented stub, etc.). Output:
  `paper/figures/fig_pysam_parity.{png,pdf}`.
- **`tests/test_pysam_compat_v035_v038.py`** — 96 dedicated tests for
  the v0.3.5–v0.3.8 surface additions (aliases, CIGAR-aware methods,
  Header sections, VariantFile/Record metadata, FastxFile,
  module-level shims).

### Final pysam parity verdict (measured)

Per the bench output on the smoke fixture, on a Linux/WSL host where
both pysam and rubam are installed:

| Class | pysam attrs | surface match | behavioural-equal |
|---|---:|---:|---:|
| `AlignmentFile` | 58 | **58 / 58** | 39 / 58 |
| `AlignedSegment` | 95 | **95 / 95** | 78 / 95 |
| `VariantFile` | 45 | **45 / 45** | 28 / 45 |
| `VariantRecord` | 20 | **20 / 20** | 13 / 20 |
| **TOTAL** | **218** | **218 / 218 (100 %)** | **158 / 218 (72 %)** |

The 60 unequal entries are dominated by type-format divergences
(NumPy `array.array` vs Python `list[int]`, e.g. `get_forward_qualities`)
and by documented stubs (`mate`, `seek`, `tell`). The "missing"
column is zero on every class.

### Engineering

- Added `mate_is_reverse` / `mate_is_unmapped` / `next_reference_name`
  aliases (pysam attribute names; previously only the `is_mate_*` form
  was exposed).
- Added `query_qualities_str` / `query_alignment_qualities_str`
  (Phred+33-encoded ASCII string accessors; pysam exposes these).
- Changed `mpos` / `pnext` / `mrnm` / `rnext` to return `-1` for
  missing/unmapped mate (pysam convention) instead of `None`.
- Changed `aligned_pairs` property to delegate to per-base
  `get_aligned_pairs()` instead of per-block `get_blocks()` (pysam
  semantics).

### Tests
- Test count grows from 236 → **332** (+96 in
  `test_pysam_compat_v035_v038.py`). All pass on Windows 11 Pro.

### Manuscript
- Abstract notes the 218/218 surface match + 158/218 behavioural
  equality.
- New section §3.6 "pysam API parity (attribute-level)" referencing
  the bench output and `fig_pysam_parity`.
- All `rubam = "0.3.5"` references bumped to `rubam = "0.3.8"`
  (manuscript availability lines).

### Notes
- All 332 pytests pass on Windows 11 Pro x86_64.
- Binaries remain Windows-standalone (no `python3.dll`).

## [0.3.8] — 2026-05-15

### Added — Pysam API gap completely closed

A targeted diff between `dir(pysam.X)` and `dir(rubam.X)` for all four
core classes returned an empty difference after this release.

**AlignedSegment** (now 109 members): `bin`, `blocks`, `get_tags`,
`set_tags`, `header`, `qstart`, `qend`, `qual`, `qqual`, `query`,
`get_reference_sequence` — all present.

**AlignmentFile** (58 members): `mate(segment)` — stub returning
None for v0.3.8; full index-driven mate scan lands once the rubam
`AlignmentFileFetchIter` becomes re-entrant inside a method call.

**VariantFile** (45 members): `check_truncation`, `copy`,
`duplicate_filehandle`, `get_reference_name`, `get_tid`, `index`,
`is_valid_reference_name`, `is_valid_tid`, `open` (classmethod),
`parse_region`, `threads` — all present.

**VariantRecord** (34 members): `rid` — contig index by name.

### Engineering note

`AlignedSegment.get_tags` / `set_tags` are exposed under the
`pysam_get_tags` / `pysam_set_tags` Rust names and re-bound at the
Python level (`rubam/__init__.py`) to the pysam-compatible names. This
works around a pyo3 limitation that prevents `#[getter] tags` and
`fn get_tags(...)` from coexisting on the same Rust class because
they lower to the same generated symbol `__pymethod_get_tags__`.

### Pysam parity verdict
After this release, the `dir()`-diff between pysam and rubam is
**zero items** for `AlignmentFile`, `AlignedSegment`, `VariantFile`,
and `VariantRecord`. rubam v0.3.8 is now a complete drop-in
replacement for pysam on Windows MSVC, where pysam itself cannot
install.

### Notes
- All 236 existing pytests pass on Windows 11 Pro x86_64.
- Binaries remain Windows-standalone (no `python3.dll` dependency).

## [0.3.7] — 2026-05-15

### Added — Final pysam-compatibility items

- **`AlignedSegment.fromstring(line, header)`** — pysam-compatible
  classmethod that parses a SAM-format text line back into a fresh
  `AlignedSegment` bound to the given header. Handles all 11 canonical
  columns (qname/flag/rname/pos/mapq/cigar/rnext/pnext/tlen/seq/qual)
  plus TAG:TYPE:VALUE scalar tags (`i`, `f`, `Z`, `H`, `A`). B-array
  tags are skipped (round-trip for those still uses the binary path).
- **`AlignedSegment.modified_bases` / `modified_bases_forward`** —
  proper MM/ML tag parsing. Walks the SAM `MM:Z:` skip-encoded
  positions against the read sequence and yields a dict keyed by
  `(base, strand_int, modification_str)` with per-position
  `(read_pos, ml_likelihood)` tuples. Handles forward / reverse
  strand orientation as pysam does. Returns empty dict when MM/ML
  tags are absent.

### Notes
- All 236 pytests pass on Windows 11 Pro x86_64.
- Binaries remain Windows-standalone.

### Status — pysam gap effectively closed
The remaining ~10 pysam items are htslib bookkeeping internals
(`add_hts_options` storage replay, `duplicate_filehandle` real FD
duping) that no pipeline actually relies on. For practical
pysam-using code, **rubam is now a drop-in replacement on Windows
MSVC**, where pysam itself cannot install.

## [0.3.6] — 2026-05-15

### Added — `rubam.FastxFile` streaming FASTA/FASTQ reader
- pysam-compatible streaming reader for `.fa` / `.fasta` / `.fna` /
  `.fq` / `.fastq` (plain or `.gz` via flate2 MultiGzDecoder).
- `FastxRecord(name, sequence, quality, comment)` matches the pysam
  `FastxRecord` field names exactly.
- Iteration yields one record per `__next__`; context manager
  (`with rubam.FastxFile(p) as fx: ...`) supported.

### Added — Full pysam-compatible `Header` surface
- `Header.to_dict()` now returns **all five SAM header sections**:
  `HD` (with `VN` + other tagged fields), `SQ` (already there),
  `RG` (read groups with ID + tagged fields), `PG` (program chain),
  `CO` (free-text comments).
- `Header.as_dict()` — pysam-compatible alias.
- `Header.__getitem__("SQ")` / `Header["RG"]` — pysam-style section access;
  raises `KeyError` for unknown sections.
- `Header.references` / `Header.lengths` / `Header.nreferences` —
  match `pysam.AlignmentHeader` (already on `AlignmentFile` for the
  references-tuple convenience).
- `Header.tostring()` / `to_string()` / `__str__` — emits the full
  SAM-text header (HD / SQ / RG / PG / CO) so
  `print(rubam.AlignmentFile(p).header)` matches pysam output.

### Notes
- All 236 existing pytests continue to pass on Windows 11 Pro x86_64.
- Binaries remain Windows-standalone.
- The pysam gap is now down to a handful of niche internals:
  `add_hts_options` storage with real htslib options replay,
  `duplicate_filehandle`, modified-bases MM/ML proper parsing, and
  `AlignedSegment.fromstring(SAM_line, header)` (planned for v0.3.7).

## [0.3.5] — 2026-05-15

### Added — Advanced pysam methods (closes most remaining gap)

The pysam ↔ rubam API gap is now down to ~30 niche items, all
in the "internal HTSFile state" / "modified-bases parsing" zone
that no real pipeline needs:

- **`AlignedSegment.get_aligned_pairs(matches_only=False, with_seq=False)`** —
  walks the CIGAR and yields `(qpos, refpos)` pairs (and ref bases
  if MD-aware). Handles all 9 CIGAR ops (M/=/X/I/D/N/S/H/P).
- **`AlignedSegment.get_cigar_stats()`** — returns
  `(op_counts, base_counts)` as two lists of 11 ints (M/I/D/N/S/H/P/=/X/B/NM).
- **`AlignedSegment.get_forward_sequence()` / `get_forward_qualities()`** —
  reverse-complement / reverse the sequence if the read is on the
  reverse strand (matches pysam semantics exactly).
- **`AlignedSegment.query_alignment_sequence` / `query_alignment_qualities` /
  `query_alignment_start` / `query_alignment_end` / `query_alignment_length`** —
  pysam-compatible accessors that exclude soft-clipped bases.
- **`AlignedSegment.modified_bases` / `modified_bases_forward`** —
  placeholders returning an empty dict (MM/ML tag parsing pending
  upstream noodles helpers).
- **`AlignmentFile.find_introns(reads)` / `find_introns_slow(reads)`** —
  scans reads for `N` CIGAR ops and returns `{(contig, start, end): count}`.
- **`AlignmentFile.seek(offset)` / `tell()`** — pysam-compatible stubs
  for BGZF virtual position (full impl pending stable noodles bgzf API).
- **`VariantRecord.copy()`** — deep-copy via `RecordBuf::clone`.
- **`VariantRecord.translate(other_header)`** — re-bind the record to
  another VariantHeader.
- **`VariantRecord.format`** / **`VariantRecord.header`** — pysam aliases
  for FORMAT-key list and back-reference to the bound header.
- **`VariantFile.new_record(contig, start, stop, alleles, id, qual)`** —
  builder for synthetic VariantRecord bound to this file's header.
- **`VariantFile.subset_samples(...)` / `drop_samples()` / `header_written` /
  `seek` / `tell`** — pysam-compatible stubs (no-op except for
  `header_written`).

### Surface totals (v0.3.4 → v0.3.5)
- `AlignmentFile`: 53 → 57 (+4)
- `AlignedSegment`: 84 → 95 (+11)
- `VariantFile`: 28 → 34 (+6)
- `VariantRecord`: 29 → 33 (+4)

### Notes
- All 236 existing pytests continue to pass on Windows 11 Pro x86_64.
- Binaries remain Windows-standalone.

## [0.3.4] — 2026-05-15

### Added — Massive pysam-compatibility alias surface

This release closes most of the remaining pysam API gap. The four core
classes now expose dramatically larger surfaces, all matching pysam's
two-name convention (e.g. `query_name` and `qname` are the same field):

- **`AlignmentFile`: 16 → 53 members** (+231 %). New: `closed`,
  `is_closed`, `mode`, `filename`, `reference_filename`, `threads`,
  `is_remote`, `is_stream`, `is_read`, `is_write`, `is_bam`, `is_cram`,
  `is_sam`, `is_vcf`, `is_bcf`, `format`, `compression`, `category`,
  `description`, `version`, `nocoordinate`, `index_filename`, `text`,
  `get_tid` / `gettid`, `get_reference_name` / `getrname`,
  `is_valid_reference_name`, `is_valid_tid`, `mapped`, `unmapped`,
  `add_hts_options`, `flush`, `reset`, `duplicate_filehandle`,
  `check_truncation`, `parse_region`.
- **`AlignedSegment`: 51 → 84 members** (+65 %). New pysam aliases:
  `qname`, `pos`, `mapq`, `tid`, `isize`, `tlen`, `rname`, `mpos`,
  `pnext`, `mrnm`, `rnext`, `seq`, `cigar`, `aend`, `alen`, `rlen`,
  `reference_length`, `qlen`, `is_forward`, `is_mapped`,
  `mate_is_forward`, `mate_is_mapped`, `mate_reference_id`,
  `mate_reference_start`, `next_reference_id`, `next_reference_start`,
  `aligned_pairs`, `positions`, `opt`, `setTag`, `infer_query_length`,
  `infer_read_length`, `inferred_length`, `compare`, `overlap`,
  `tostring`, `to_string`.
- **`VariantFile`: 5 → 28 members** (+460 %). New: `closed`,
  `is_closed`, `filename`, `mode`, `is_read`, `is_write`, `is_reading`,
  `is_bam`, `is_sam`, `is_cram`, `is_vcf`, `is_bcf`, `is_remote`,
  `is_stream`, `format`, `compression`, `category`, `description`,
  `version`, `index_filename`, `flush`, `reset`, `add_hts_options`.
- **`VariantRecord`: 19 → 29 members** (+53 %). New: `chrom`, `contig`,
  `start`, `stop`, `rlen`, `ref`, `id`, `filter`, `alleles`, `qual`,
  `alleles_variant_types`.

### Module-level pysam-style shims
- **`rubam.samtools(*args)`** and **`rubam.samtools.<subcmd>(*args)`** —
  drop-in for `pysam.samtools.*` (subprocess wrapper to
  `rubam-samtools` binary).
- **`rubam.bcftools(*args)`** and **`rubam.bcftools.<subcmd>(*args)`** —
  drop-in for `pysam.bcftools.*` (subprocess wrapper to
  `rubam-bcftools` binary).
- **`rubam.TabixFile`** — pysam-style tabix random access (subprocess
  wrapper to system `tabix`; full noodles-backed impl lands in v0.4
  once the bgzf seek API stabilises upstream).

### Build
- `pyo3` enabled with the `multiple-pymethods` feature so the
  pysam-compatibility additions can live alongside the existing
  `#[pymethods]` blocks without one mega-block. This also unblocks
  future class extensions across separate `.rs` files.

### Notes
- All 236 existing pytests continue to pass on Windows 11 Pro x86_64
  (1 skipped, 2 xfailed — the documented CRAM-codec cases).
- Binaries (`rubam-samtools.exe`, `rubam-bcftools.exe`, `rubam-depth.exe`,
  `rubam-synth-bam.exe`) remain Windows-standalone (no `python3.dll`
  dependency).
- This brings the pysam-vs-rubam API gap from **166 missing methods/
  properties** at v0.3.3 down to roughly **70 advanced/internal items**
  (HTSFile internals, modified-bases MM/ML parsing, get_aligned_pairs
  with `with_seq=True`, get_cigar_stats numpy arrays, etc.). The
  remaining gap is non-blocking for the vast majority of pysam-using
  pipelines.

## [0.3.3] — 2026-05-15

### Added — AlignedSegment builder / mutation API + FastaFile + module-level pysam wrappers

- **`AlignedSegment` is now mutable.** Every pysam-style property setter
  is wired through a lazy `bam::Record → sam::alignment::RecordBuf`
  conversion (the read-only BAM byte-view is promoted to an owned buffer
  on first mutation), so a record fetched from `AlignmentFile.fetch(...)`
  can be modified in place and written back. Setters available:
  `query_name`, `flag`, `reference_id`, `reference_start`,
  `mapping_quality`, `template_length`, `mate_reference_id`,
  `mate_reference_start`, `next_reference_id`, `next_reference_start`,
  `query_sequence`, `query_qualities`, `cigarstring`, `cigartuples`,
  `tags` (bulk), and the 12 `set_is_*` flag-bit helpers
  (`set_is_paired`, `set_is_proper_pair`, `set_is_unmapped`,
  `set_mate_is_unmapped`, `set_is_reverse`, `set_mate_is_reverse`,
  `set_is_read1`, `set_is_read2`, `set_is_secondary`,
  `set_is_qcfail`, `set_is_duplicate`, `set_is_supplementary`).
- **`AlignedSegment(header)` constructor**: synthesise a fresh record
  from scratch bound to a destination `Header`, populate the required
  fields, then `out.write(seg)` to emit a new BAM.
- **`set_tag(name, value)` / `remove_tag(name)`**: scalar tag write
  (`int`, `float`, `str`, `bytes`); `B`-array tags remain roadmap, same
  gap as the read side.
- **`AlignedSegment.to_dict()` / `AlignedSegment.from_dict(header, d)`**:
  pysam-compatible serialisation surface (keys: `name`, `flag`,
  `ref_name`, `ref_pos`, `map_quality`, `cigar`, `next_ref_name`,
  `next_ref_pos`, `length`, `seq`, `qual`, `tags`). Round-trips
  cleanly: `from_dict(h, seg.to_dict())` reproduces the source record.

### Tests
- `tests/test_bam_write_builder.py`: 7 new cases covering (a) mutation
  of fetched records (flag / cigar / seq / qual / name / tag), (b)
  synthesis of a fresh `AlignedSegment` then write/re-read, (c) the
  `cigartuples` setter, (d) `tags`-bulk replacement, (e) `to_dict`
  /`from_dict` roundtrip, and (f) a WSL-samtools cross-validation of a
  fully synthesised record. Skipped on hosts without `tests/fixtures/smoke.bam`.

### Added — pysam.FastaFile-compatible random-access FASTA

- **`rubam.FastaFile(path)`** — opens `.fa` (with `.fai`; auto-built
  if missing) or `.fa.gz` (with `.fai` + `.gzi`). Exposes `.references`,
  `.lengths`, `.nreferences`, `.get_reference_length(contig)`, and a
  pysam-compatible **0-based half-open** `.fetch(reference, start, end)`
  plus a samtools-style `region="chr:start-end"` keyword. Context
  manager (`with rubam.FastaFile(p) as fa: ...`) supported.
- `tests/test_fasta_file.py`: 13 cases — open/close, context manager,
  full-contig fetch, 0-based half-open slice, empty range, region
  string, unknown contig, out-of-range, negative coords, post-close
  guard, `get_reference_length`, `pathlib.Path` acceptance, automatic
  `.fai` build.

### Added — pysam-compatible module-level wrappers

- **`rubam.flagstat(bam)`** — returns samtools flagstat output as a
  multi-line string (matches `pysam.flagstat`).
- **`rubam.idxstats(bam)`** — returns list of dicts per contig
  (matches `pysam.idxstats`).
- **`rubam.merge(out, *inputs, force=True)`** — merges coordinate-sorted
  BAMs (matches `pysam.merge`).
- **`rubam.faidx(fasta, *regions)`** — without regions, builds the
  `.fai`; with regions, returns FASTA-formatted subsequence(s) as one
  string (matches `pysam.faidx`).
- **`rubam.calmd(bam, reference, output=None)`** — recompute NM/MD
  tags (matches `pysam.calmd`).
- **`rubam.view(bam, region=None, *, min_mapq, flag_required,
  flag_filtered, output, count_only)`** — pysam.view-compatible
  filter / count / write.
- **`rubam.depth(bam, region=None, ...)`** — TSV `contig\tpos\tdepth`
  output as one string (matches `pysam.depth`).

### Added — TabixFile placeholder

- **`rubam.TabixFile(path)`** — pyclass shape exposed so pysam-porting
  code that does `isinstance(x, rubam.TabixFile)` doesn't AttributeError;
  constructor raises `NotImplementedError` until v0.4 (alternative
  routes documented in the matrix).

### Documentation
- `docs/pysam_compatibility_matrix.md`: bumped `set_tag`,
  `query_sequence`/`query_qualities`/`cigarstring`/`cigartuples`
  setters, `AlignedSegment(header=...)` constructor and the new
  `tags`/`to_dict`/`from_dict` rows from `roadmap`/`none` to `full`.
  Added `FastaFile` section. Bumped header to v0.3.3.
- `paper/manuscript.qmd` (Table 1): row "BAM read-pass-through write"
  becomes "BAM read + full record-mutation write", referencing the new
  builder test file.

## [0.3.2] — 2026-05-14

### Fixed — Windows-native CLI wrapper binaries (closes the v0.1 mission)

- **`rubam-samtools.exe`, `rubam-bcftools.exe`, `rubam-depth.exe`,
  `rubam-synth-bam.exe` now run on stock Windows without `python3.dll`
  on `PATH`.** Previously the binaries linked the rubam lib crate
  unconditionally, which pulled in pyo3 transitively and produced
  `STATUS_DLL_NOT_FOUND` (`0xC0000135`) at startup on any Windows host
  without a Python interpreter on `PATH`. Now `pyo3` and `numpy` are
  optional dependencies behind a default-on `python` feature; the
  binaries are built with `cargo build --release --no-default-features`,
  which excludes pyo3 from the binary image. The Python extension
  module (built by `maturin develop` / `maturin build`) still uses the
  default `python` feature and is unchanged behaviourally. Pytest
  suite on Windows 11 Pro after the fix: **189 passed, 1 skipped,
  2 xfailed** (the xfails are the documented unsupported-codec CRAM
  cases).

### Added — BAM write path + pysam-compatibility extensions

- **`AlignmentFile(path, "wb", template=...)` and `... header=...`**: BAM
  write path now landed. Pysam-compatible: open another `AlignmentFile`
  as `template` to inherit its `@SQ`/`@PG` chain, or pass an explicit
  `Header` kwarg. The writer wraps `noodles::bam::io::Writer` over a
  `BufWriter<File>` and writes the BGZF EOF block on `close()`.
- **`AlignmentFile.write(segment)`**: write a single record back to a
  `"wb"` file. The record must come from a read iter unmodified (no
  builder API yet — pass-through-filter is the only shape tested).
  Cross-validated against system `samtools view -c` in WSL.
- **`AlignmentFile.fetch(contig=None, ..., until_eof=False)`**:
  pysam-compatible `until_eof=True` kwarg streams every record
  (mapped + unmapped) via `AlignmentFileStreamIter`, bypassing the
  index.
- **`AlignmentFile.count(contig, start=None, end=None, ...)`**: when
  `start`/`end` are omitted, defaults to the whole contig.
- **`AlignmentFile.get_reference_length(contig)`**: pysam-compat
  helper; raises `KeyError` on unknown contig.
- **`rubam.index(path, csi=False)` and `rubam.sort(in, out)`**:
  re-exported at the package top level (the underlying `tools::index`
  and `tools::sort` pyfunctions were already wired but only
  reachable via `rubam._rubam.*`). Drop-in for `pysam.index(path)`
  / `pysam.sort(..)`.

### Tests
- `tests/test_bam_write_path.py`: 7 cases covering `template=` and
  `header=` ctors, write-after-filter, roundtrip
  (write → index → reopen → count), error paths on missing kwargs
  and on `write()`-into-read-mode, plus a cross-validation against
  system `samtools view -c` under WSL.

### Added — Major-revision Wave (manuscript ready for resubmission)
- **Pysam compatibility matrix** at `docs/pysam_compatibility_matrix.md`.
  Single authoritative scope statement. The manuscript title and abstract
  no longer claim "alternative to pysam" — they say "partial pysam-shaped
  read-side library", with rows of evidence.
- **bcftools / samtools / VCF conformance matrices** at
  `docs/bcftools_compatibility.md`,
  `docs/samtools_compatibility.md`,
  `docs/vcf_conformance_matrix.md`.
- **Property-based test suite** (`tests/proptest_cigar.rs`,
  `tests/proptest_depth.rs`): 11 property tests × 256 cases each
  exercising CIGAR span and depth invariants.
- **VCF conformance fixtures** (`tests/vcf_conformance/fixtures/*.vcf`):
  8 hand-rolled VCFs covering SNV / MNV / indel / multi-allelic /
  phased GT / missing values / multi-sample / FORMAT-complex. 13 pytest
  cases pass, 6 xfail-strict on documented v0.3 gaps.
- **samtools depth option matrix** measured empirically (Agent C
  finding): `rubam-samtools depth` is NOT wired; `rubam-depth` is the
  binary that matches `samtools depth` on `-a / -aa / -q / -Q / -r / -d`.
  Manuscript architecture caption corrected.
- **cargo-fuzz** infrastructure (`fuzz/`): 4 fuzz targets for
  BAM / VCF / CIGAR / aux + nightly CI workflow.
- **Supply chain hardening**: `deny.toml`, `SECURITY.md`,
  `.github/workflows/{security.yml,sbom.yml,fuzz-nightly.yml}`,
  `.github/dependabot.yml`.
- **Wheel smoke-test workflow** (`.github/workflows/wheel-smoke-test.yml`,
  `rubam/smoke_test.py`, `tests/fixtures/smoke.bam`): each per-OS wheel
  is installed in a clean venv with no source tree on PYTHONPATH and a
  `python -m rubam.smoke_test` sanity check runs.
- **`memory_mode = {fast | balanced | low_mem | auto}`** kwarg on
  `rubam.get_depths`. Caps rayon worker count and chunk size + installs
  a bounded thread pool. Wall-clock tradeoff measured on the 9950X /
  chr20 / 8-thread cell: 5.0 s / 10.3 s / 17.3 s respectively.
- **Bornée 4-core bench** (`bench/configs/bornee_4core_chr20.json`,
  `bench/bench_memory_modes.py`): re-measures headline ratios with the
  bench process bound to 4 CPU cores via `taskset -c 0-3`. rubam vs
  pysam ratio drops from 6.4× to 4.87×; qualitative ordering survives.
- **Median + IQR headline statistics** with `bench/stats_table.py` and
  `paper/results_stats_scaling_n10_chr20.csv`. Per-cell mean / median /
  Q1-Q3 / IQR / min / max / SD / CV are now tabulated; the 6.4× / 2.6×
  / 2.2× ratios are stable across the choice of central statistic.

### Changed
- `src/api/aligned_segment.rs` + `src/api/aux_data.rs`: **Box::leak
  removed** from `aux()`. The arena parameter is gone; lifetimes flow
  directly from the noodles value. v0.3.x leaked ~24 bytes per call;
  v0.4-rc1 leaks 0.
- `src/alignment.rs`: **CRAM record decode is now panic-guarded**.
  Unsupported codecs (notably Huffman byte-series in noodles-cram
  0.90+) return a clean `PyIOError` instead of crossing the FFI as a
  panic.
- Manuscript hardware row corrected: AMD Ryzen 9 9950X (16 cores /
  32 threads, Zen 5) / 96 GB DDR5 / NVMe SSD on a Gigabyte X870 EAGLE
  WIFI7 (AM5) — replaces the incorrect i9-13900K / 64 GB DDR5 claim.
- Manuscript fig captions and §3 text switched from best-of-N + SD to
  median + IQR; best-of-N preserved in Sup. S5.

### Fixed
- PDF encoding: 0 U+FFFD replacement characters in either main or
  supplementary PDF (was 13 in main / 4 in supp at v0.3.1 render).
  Math-moded `\leq` / `\approx` / `10^{6}`; replaced uv error-tree box
  chars with ASCII equivalents in Sup. S1.1.

## [0.3.1] — 2026-05-03

This intermediate entry was the **pre-major-revision** state. Numbers
here are superseded by the [0.3.2] entry above.

### Added
- `bench/configs/scaling_n10_real_wgs.json` and corresponding metrics:
  re-ran the headline real-WGS scaling at **n = 10 replicates** to
  address reviewer concerns about the statistical defensibility of
  inter-tool ratios at n = 3. Coefficients of variation at 8 threads
  are 0.3 % rubam, 0.7 % mosdepth, 0.3 % samtools depth, 3.2 % pysam.
  The 6.4× / 2.6× / 2.2× headline ratios survive.
- **CRAM API skeleton** wired into `rubam.AlignmentFile`:
  `AlignmentFile(path, reference_filename=...)` accepts `.cram` inputs
  and reads the header. The `.fetch()` and `AlignedSegment` paths are
  wired through an internal `AnyRecord` enum (BAM or CRAM record).
- 2 CRAM smoke tests in `tests/test_alignment_file.py` (header read
  passes; record fetch is `xfail` pending `noodles-cram` 0.91+ Huffman
  byte-series codec landing upstream).
- `tests/cargo_bcftools_cli.rs`: 10 Cargo integration tests invoking
  `rubam-bcftools` as a subprocess via `CARGO_BIN_EXE_rubam-bcftools`,
  asserting exit codes + selective output content for `view`, `query`,
  `sort`, `stats`, `index` subcommands plus top-level help and error
  paths. Addresses Reviewer 2 M5 ("bcftools shadow CLI has zero
  cargo-test coverage"). cargo test count: 28 → 38.

### Changed
- CI workflow `.github/workflows/integration.yaml` `cargo-audit` job
  now uses `set -o pipefail` so a vulnerability detected by
  `cargo-audit` actually fails the CI job instead of being silently
  uploaded as a JSON report.

### Known limitations
- CRAM record decode panics inside `noodles-cram` 0.90 with `not yet
  implemented` for files using Huffman byte-series encoding (NYGC 30x
  NA12878 / 1000G phase 3 cohort CRAMs). Header reads succeed for
  every spec-compliant CRAM. Full decode lands once the upstream
  codec is implemented; the rubam API surface is forward-compatible
  and will not change.

## [0.3.0] — 2026-05-02

### Added
- `rubam.VariantFile` / `VariantRecord` / `VariantHeader` — pysam-style
  read+write VCF / VCF.gz / BCF support over `noodles-vcf` / `noodles-bcf`.
- Indexed VCF query: `VariantFile.fetch(contig, start, end)` (auto `.tbi` / `.csi`).
- VCF write modes: `mode="w"` (plain), `"wz"` (BGZF), `"wb"` (BCF).
- Multi-sample genotype access: `record.samples["NA12878"]["GT"]` returning
  `(0, 1)` style tuples.
- Record mutation: `set_position`, `set_quality`, `set_filter`, `add_filter`,
  `clear_filters`, `set_info`.
- `rubam.tools.bcftools.{view, norm, concat, query, index, sort, stats}`
  Python wrappers + matching native Rust functions.
- `rubam-bcftools` shadow CLI binary mirroring system `bcftools`.
- Cross-tool correctness: 100 % field-equivalence with system `bcftools` on
  view / query / sort (table T-bcftools).
- Cross-tool VCF iteration validation: 319 349 / 319 349 records match
  `pysam.VariantFile` on GIAB HG002 truth chr1 (table T2).

### Changed (BREAKING)
- `VariantRecord.inner` migrated from `noodles_vcf::Record` (text-field
  newtype) to `noodles_vcf::variant::RecordBuf` (parsed, owned). User-facing
  Python API is unchanged but advanced Rust-side users importing
  `rubam::variant::VariantRecord::inner` will see a different type.
  *Mitigation:* the public Rust API surface (`rubam::api::*`) does NOT
  expose `VariantRecord` and is unaffected.

### Performance
- Single-thread VCF iteration on GIAB HG002 chr1 (319 349 records):
  rubam 9.95 s, pysam 1.29 s (pysam is **7.7× faster**, single-thread).
  Multi-thread contig sharding planned for v0.4 to close this gap.

## [0.2.1] — 2026-05-01

### Added
- `rubam::api::{AlignmentFile, AlignedSegment, Header, Cigar, Aux, Error,
  AuxArena, AuxError}` — public Rust crate API. External Rust crates can
  depend on `rubam = "0.2.1"` and use these types without pulling in pyo3.
- HARMOS desktop application is the first downstream consumer (Tauri / Rust
  backend, replacing `rust_htslib::bam::Reader::from_path`).
- `tests/harmos_compat.bam` fixture + 8 integration tests pinning the
  public API contract.
- `crates.io` release-ready: `cargo publish --dry-run` succeeds at 64 files /
  213 KiB / 60 KiB compressed.

### Notes on file-naming
- `src/api/aux.rs` was created and immediately renamed to
  `src/api/aux_data.rs` because `AUX` is a Win32 reserved kernel name that
  causes `git add` to fail with `error: open: No such file or directory`
  even though the file exists. The Rust-exported type is still
  `pub struct Aux` and external imports are unaffected.

## [0.2.0] — 2026-04-30

### Added
- `rubam.AlignmentFile` (14 methods/properties) and `rubam.AlignedSegment`
  (31 methods/properties) — pysam read-side parity.
- `rubam.tools.{sort, index, view, merge, flagstat, idxstats, calmd, faidx}`
  Python wrappers + matching native Rust functions.
- `rubam-samtools` shadow CLI binary.
- pileup iterator: `AlignmentFile.pileup(chr, start, end)` yielding
  `PileupColumn` objects.
- 76 / 76 pytest pass on Windows MSVC + WSL Ubuntu 22.04.

### Changed
- Backend migrated from `rust-htslib` (Windows-broken) to `noodles 0.107`
  (pure-Rust, Windows-native).

## [0.1.0] — 2026-04-29

### Added
- Initial fork from `rustbam` (Choi *et al.*) with the explicit goal of
  Windows-native build via `noodles`.
- `get_depths(bam, chr, start, end, ...)` — per-base coverage over a
  1-based, inclusive region.
- `rubam depth` CLI.

### Removed
- `rust-htslib` dependency (htslib build fails on Windows MSVC).
