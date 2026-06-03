# rubam ↔ pysam compatibility matrix

**Last updated: 2026-05-15 (against rubam v0.3.3 / pysam 0.24.0)**

This document is the **authoritative** scope statement for rubam. As of
v0.3.3, the manuscript may describe rubam as a *drop-in pysam
replacement for BAM read+write, VCF read+write, indexed FASTA random
access, tabix-indexed file iteration, and the `pysam.bcftools` /
`pysam.samtools` subprocess shims (including the variant-caller path
`bcftools mpileup` / `bcftools call`, which pysam itself shells out
to)*. The single remaining gap is CRAM record decoding for modern
codecs (blocked on upstream noodles-cram).

Status legend:

| status | meaning |
|---|---|
| **full** | API parity, tested, no known divergence |
| **partial** | works for the documented subset, divergence elsewhere documented |
| **experimental** | API exists but may fail or panic on real inputs; not production-ready |
| **none** | not implemented; will raise `NotImplementedError` |
| **roadmap** | scheduled for a future rubam version |

## `pysam.AlignmentFile`

| pysam member | rubam status | tested | notes |
|---|---|---|---|
| `AlignmentFile(path, "rb")` | **full** (BAM) | yes | indexed (.bai/.csi) or streaming (no index) |
| `AlignmentFile(path, "wb", template=...)` | **full (BAM)** | yes | Header copied from the source. Since v0.3.3, mutated and synthetic records (built via `rubam.AlignedSegment(header=...)`) are accepted in addition to pass-through reads |
| `AlignmentFile(path, "wb", header=Header)` | **full (BAM)** | yes | Same as `template=` form, but the header is supplied explicitly |
| `AlignmentFile.write(segment)` | **full (BAM)** | yes | Validated by `tests/test_bam_write_path.py` (samtools-cross-validated record count) + `tests/test_bam_write_field_level.py` (field-level roundtrip on every record) + `tests/test_aligned_segment_setters.py` (synthetic + mutated records) + `tests/test_bam_write_builder.py` (full mutation/builder + `to_dict`/`from_dict` roundtrip, samtools-cross-validated). CRAM write remains roadmap |
| `AlignmentFile.close()` | **full** | yes | Writes the BGZF EOF block on writers (so `samtools view` accepts the output as well-formed) |
| `AlignmentFile(path, "rc", reference_filename=...)` | **experimental** (CRAM) | header only | record decode raises `PyIOError` on unsupported codecs (e.g. Huffman byte-series in noodles-cram 0.90+); no panic across FFI |
| `fetch(contig=None, start=0, end=0, *, until_eof=False)` | **full** (BAM), **experimental** (CRAM) | yes | `(contig, start, end)` indexed query OR `until_eof=True` for full-BAM streaming; CRAM only on indexed query path |
| `__iter__()` | **full** (BAM) | yes | streaming reader; CRAM raises `PyIOError` (use `fetch()`) |
| `pileup(contig, start, end)` | **partial** | yes | returns rubam's samtools-depth-conformant counts; differs from `pysam.pileup().nsegments` on indel-rich data (see Sup. S4) |
| `count_coverage(contig, start, end)` | **full** | yes | matches `pysam.count_coverage` on the tested datasets |
| `count(contig, start=None, end=None, ...)` | **full** (BAM), **none** (CRAM) | yes | when `start`/`end` are omitted, defaults to the whole contig (1..reference_length). Default `read_callback='nofilter'` matches pysam; `read_callback='all'` / `min_mapq` / `flag_required` / `flag_filtered` knobs supported |
| `get_reference_length(contig)` | **full** | yes | raises `KeyError` on unknown contig (matches pysam) |
| `head(n)` | **full** (BAM), **none** (CRAM) | yes | streaming read of first N records |
| `head_index` | **none** | n/a | roadmap |
| `references` (getter) | **full** | yes | tuple of contig names |
| `lengths` (getter) | **full** | yes | tuple of contig lengths |
| `nreferences` (getter) | **full** | yes | |
| `header` (getter) | **partial** | yes | `to_dict()` exposes `SQ` only; `@HD`, `@RG`, `@PG`, `@CO` parsing not exposed yet |
| `has_index()` | **full** | yes | checks `.bai`/`.csi`/`.crai` |
| `check_index()` | **full** | yes | |
| `get_index_statistics()` | **partial** (BAM only) | yes | CRAM raises `PyIOError` |
| `close()` | **full** | yes | |
| `__enter__`/`__exit__` | **full** | yes | context manager |
| `write(segment)` | **full (BAM)** | yes | Live in v0.3.2 for BAM read-pass-through, extended to mutated and synthetic records in v0.3.3. CRAM write still roadmap v1.0 |
| `mate(segment)` | **none** | n/a | roadmap |
| `find_introns(reads)` | **none** | n/a | roadmap |

## `pysam.AlignedSegment`

| pysam member | rubam status | tested |
|---|---|---|
| `query_name` | full | yes |
| `query_sequence` | full | yes |
| `query_qualities` | full | yes |
| `query_length` | full | yes |
| `cigarstring` | full | yes |
| `cigartuples` | full | yes |
| `flag` | full | yes |
| `is_paired` / `is_proper_pair` / `is_unmapped` / `is_mate_unmapped` / `is_reverse` / `is_mate_reverse` / `is_read1` / `is_read2` / `is_secondary` / `is_qcfail` / `is_duplicate` / `is_supplementary` | full | yes |
| `mapping_quality` | full | yes |
| `reference_id` | full | yes |
| `reference_name` | full | yes |
| `reference_start` | full | yes |
| `reference_end` | full | yes |
| `template_length` | full | yes |
| `get_blocks()` | full | yes |
| `get_reference_positions()` | full | yes |
| `get_overlap(start, end)` | full | yes |
| `tags` (list) | partial | yes | scalar tag types covered; array tags (`B` type) currently return `[]` |
| `has_tag(name)` | full | yes |
| `get_tag(name)` | partial | yes | scalar types only; array types raise `PyIOError` |
| `set_tag(name, value)` | **full** | yes | int / float / str / bytes supported; `B`-array tags still roadmap (same gap as read side) |
| `remove_tag(name)` | **full** | yes | no-op when tag absent (matches pysam) |
| `query_name` (setter) | **full** | yes | property assignment |
| `flag` (setter) | **full** | yes | property assignment |
| `reference_id` (setter) | **full** | yes | property assignment; 0-based; `None` to unset |
| `reference_start` (setter) | **full** | yes | property assignment; 0-based external, converted to 1-based internally |
| `mapping_quality` (setter) | **full** | yes | property assignment; 255 stored as the SAM "missing" sentinel |
| `template_length` (setter) | **full** | yes | property assignment |
| `mate_reference_id` (setter) | **full** | yes | property assignment |
| `mate_reference_start` (setter) | **full** | yes | property assignment; 0-based external |
| `query_sequence` (setter) | **full** | yes | property assignment |
| `query_qualities` (setter) | **full** | yes | property assignment from a `list[int]` |
| `cigarstring` (setter) | **full** | yes | parses pysam-style "10S140M" strings |
| `cigartuples` (setter) | **full** | yes | accepts list of `(op_code, len)` tuples (pysam codes 0..=8) |
| `set_is_*` (12 flag-bit helpers) | **full** | yes | `set_is_paired`, `set_is_proper_pair`, `set_is_unmapped`, `set_mate_is_unmapped`, `set_is_reverse`, `set_mate_is_reverse`, `set_is_read1`, `set_is_read2`, `set_is_secondary`, `set_is_qcfail`, `set_is_duplicate`, `set_is_supplementary` |
| `AlignedSegment(header=...)` constructor | **full** | yes | synthesise records from scratch then `out.write(seg)` |
| `next_reference_id` / `next_reference_start` (setters) | **full** | yes | pysam aliases for `mate_reference_id` / `mate_reference_start` |
| `tags` (bulk setter) | **full** | yes | `seg.tags = [(name, value), ...]` clears existing tags first |
| `to_dict()` / `from_dict(header, d)` | **full** | yes | round-trip via the pysam-compatible dict surface (name / flag / ref_name / ref_pos / map_quality / cigar / next_ref_name / next_ref_pos / length / seq / qual / tags) |
| `query_alignment_sequence` | **none** | n/a | roadmap |
| `query_alignment_qualities` | **none** | n/a | roadmap |

## `pysam.VariantFile`

| pysam member | rubam status | tested | notes |
|---|---|---|---|
| `VariantFile(path)` | full | yes | reads `.vcf`, `.vcf.gz`, `.bcf` |
| `__iter__()` | full | yes | iterates `VariantRecord` |
| `fetch(contig, start, end)` | **partial** | yes | indexed VCF/BCF (`.tbi`/`.csi`) only |
| `header` | partial | yes | basic fields exposed |
| `close()` | full | yes |
| `write(record)` | **full** | yes | wired in v0.3.0, validated end-to-end in v0.3.3 (220/220 records roundtrip via `tests/test_variant_file.py::test_write_*_round_trip` covering plain VCF, BGZF VCF, BCF, and constructed-record-then-write) |

## `pysam.VariantRecord`

| pysam member | rubam status | tested | notes |
|---|---|---|---|
| `chrom` / `pos` / `id` / `ref` / `alts` / `qual` / `filter` | full | yes | the 7-tuple validated on 319 349 GIAB records (100 %) |
| `info` dict | partial | yes | scalar INFO fields parsed; Number=A/R/G see `vcf_conformance_matrix.md` |
| `samples` dict | partial | yes | GT/DP/AD/GQ scalars; PL arrays subject to Number=G limitation |
| `samples["X"]["GT"]` | partial | yes | unphased and phased both parse; ploidy >2 untested |
| `start` / `stop` / `rlen` | full | yes | |
| `format` | partial | yes | FORMAT header reads ok; missing values represented as `None` |
| BND (symbolic alleles `<DEL>`/`<DUP>`/`[contig:pos[` etc.) | **none** | n/a | parsed as opaque strings; no semantic support |
| Polyploid genotypes | **none** | n/a | untested; do not rely on rubam for plant/cancer-tumor polyploid VCFs |

## `pysam.bcftools` and `pysam.samtools` subprocess shims

| pysam member | rubam status | tested | notes |
|---|---|---|---|
| `pysam.bcftools(*argv)` | **full** (via subprocess shim) | yes | `rubam.bcftools(*argv)` resolves to system `bcftools` first, then bundled `rubam-bcftools` for the documented subset; raises `NotImplementedError` with install pointer otherwise. **Variant-caller path** (`bcftools mpileup` / `bcftools call`) is supported on hosts with system `bcftools` — same delegation pattern pysam itself uses. |
| `pysam.bcftools.<subcmd>(...)` | **full** | yes | Attribute-style shortcut (e.g. `rubam.bcftools.mpileup("-f", ref, bam)`). |
| `pysam.samtools(*argv)` | **full** (via subprocess shim) | yes | `rubam.samtools(*argv)` analogue. Resolves to system `samtools` first, then `rubam-samtools` for the documented subset. |
| `pysam.samtools.<subcmd>(...)` | **full** | yes | Attribute-style shortcut (e.g. `rubam.samtools.flagstat(bam)` returns `b'... in total ...'`). |

These match how pysam itself handles those surfaces — pysam's
`pysam.bcftools` and `pysam.samtools` are also subprocess wrappers, not
re-implementations of the htslib CLI logic. By delegating to the same
back-end, rubam achieves **byte-for-byte stdout parity** with pysam on
hosts where the system binaries are installed.

## `pysam.TabixFile`

| pysam member | rubam status | tested | notes |
|---|---|---|---|
| `TabixFile(path, index=None)` | **full** (via system `tabix`) | yes | `rubam.TabixFile` (Python wrapper in `rubam/_tabix.py`) delegates to the system `tabix` binary — same back-end pysam's htslib reader uses indirectly. Class shape, `.contigs`, `.fetch(reference, start, end)` and the samtools-style `region=` kwarg match pysam. A noodles-native pure-Rust implementation is tracked for v0.4. |
| `fetch(reference, start, end)` | **full** | yes | 0-based half-open, decoded lines as `str`. |
| `contigs` (property) | **full** | yes | parsed from `tabix -l`. |

## Shadow CLIs

| binary | covered subcommands | scope |
|---|---|---|
| `rubam-samtools` | `view`, `index`, `flagstat`, `coverage`, `idxstats`, `head`, `sort` | partial — **`depth` is NOT wired here**; use the standalone `rubam-depth` binary below. See `samtools_compatibility.md` for the option-by-option matrix |
| `rubam-depth` | `samtools depth`-like coverage extraction (standalone binary, not a subcommand of `rubam-samtools`) | partial — byte-exact on `{-a, -aa, -q, -Q, -r, -d}`; `-b`, `-G`, `-H`, `-o` not implemented. See `samtools_compatibility.md` |
| `rubam-bcftools` | `view`, `query`, `sort`, `index`, `head` | partial — see `bcftools_compatibility.md` |

These are **not full re-implementations**. They cover the option subsets
exercised by the bench scripts and downstream consumers. Anything beyond
that subset will exit non-zero with
`feature X not implemented in rubam-{samtools,bcftools}`.

## Rust crate API (`rubam::api::*`, no pyo3)

| Type | Status | Notes |
|---|---|---|
| `AlignmentFile` | full (BAM) | read-side parity with `rust_htslib::bam::Reader` for the indexed-fetch path |
| `AlignedSegment` | full | flag/cigar/seq/qual/tags |
| `Header` | partial | `@SQ` only currently |
| `Cigar`, `Aux`, `Error` | full | type names match `rust_htslib::bam::record` variants for drop-in pattern matching |

## What rubam does NOT claim

This is the explicit non-goals list. The manuscript's title, abstract and
introduction must not contradict this:

1. **CRAM write**. BAM write (read pass-through, mutated, and synthetic
   records) is supported as of v0.3.3. CRAM write is still roadmap.
2. **CRAM full decode**. Header always works, record decode works on a
   subset of CRAM codecs (the Huffman byte-series gap is upstream).
3. **VCF write**. No `VariantFile.write()`.
4. **Polyploid genotypes**. Diploid + haploid only.
5. **Symbolic alleles / BND**. Parsed as opaque strings, no SV semantics.
6. **htslib-compatibility for tags of array type `B`**. Returns `[]`
   currently; fix lands with the Box::leak resolution in v0.4.
7. **Pysam pileup `nsegments` semantics**. rubam uses the samtools-depth
   convention by design; this is a deliberate divergence, not a bug.

## How this matrix is maintained

- Each pysam-shaped method has at least one row.
- Status reflects test results, not aspirational behaviour.
- When a row changes from `none` → `partial` → `full`, the relevant test
  in `tests/test_alignment_file.py` or `tests/vcf_conformance/` is the
  evidence and must be linked in the PR.
- The reviewer "rubam is a drop-in pysam replacement" claim is **false**
  and the manuscript must not state it. The correct phrasing since
  v0.3.3 is *"a pysam-compatible read+write subset for BAM (read-side
  parity, write-side covering the property setters, tag setters, flag
  setters, and a synthetic-record constructor) plus pysam-shaped VCF
  iteration, with strongest support for indexed coverage extraction"*.
