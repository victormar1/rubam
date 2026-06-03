# rubam-bcftools ↔ system bcftools compatibility

**Last updated: 2026-05-14 (against rubam v0.3.2 / bcftools 1.x)**

This matrix is the authoritative scope of the `rubam-bcftools` CLI wrapper.
The manuscript must not say "full bcftools surface mirrored" because that
is false. The accurate phrasing is *"a subset of bcftools subcommands
covering the most common BAM/VCF inspection paths used by the benchmark
scripts and downstream consumers"*.

Status legend:

| status | meaning |
|---|---|
| **full** | rubam-bcftools matches system bcftools on the tested option set |
| **partial** | the subcommand exists but only a documented option subset works |
| **partial+caveat** | works on most inputs, with a known semantic gap |
| **none** | not implemented; the binary exits 2 with "subcommand not implemented" |
| **roadmap** | scheduled for a future release |

## Subcommand status

| bcftools subcommand | rubam status | tested records | output match | notes |
|---|---|---|---|---|
| `view` | partial | 1 279 | **100 %** field-equivalent | region filter `chr:start-end` works; `-O b` (BCF output) supported; `-O v` (VCF text) supported; sample filters (`-s`, `-S`) untested |
| `query` | partial | 134 | **100 %** | format strings: `%CHROM`, `%POS`, `%REF`, `%ALT`, `%QUAL`, `%FILTER`, `%INFO/<TAG>` for scalar tags. Array tags (Number=A/R/G) and multi-sample format strings untested |
| `sort` | partial | 1 000 | **100 %** | in-memory sort by chrom/pos; temp-file external sort for huge files is roadmap |
| `index` | partial | n/a | byte-identical CSI | builds `.csi` index for bgzipped VCF; `.tbi` (tabix) emit untested |
| `head` | full | n/a | byte-identical | first N records, mirrors system bcftools 1.18+ behaviour |
| `norm` | partial+caveat | n/a | **diverges** on Number=G | left-align + split multi-allelics works on Number=A/R; the Number=G genotype-likelihood reshape uses a simplified rule and may differ from system bcftools — open issue, not yet a regression test |
| `merge` | **none** | n/a | n/a | roadmap v0.5 |
| `concat` | **none** | n/a | n/a | roadmap v0.5 |
| `annotate` | **none** | n/a | n/a | not planned in v0.x; out-of-scope |
| `stats` | **none** | n/a | n/a | roadmap v0.5 |
| `isec` | **none** | n/a | n/a | not planned in v0.x |
| `filter` | **none** | n/a | n/a | roadmap v0.5 |
| `gtcheck` | **none** | n/a | n/a | not planned |
| `mpileup` | **none** | n/a | n/a | covered by `rubam-samtools pileup` |
| `call` | **none** | n/a | n/a | not planned (out of scope; rubam is not a variant caller) |
| `roh` | **none** | n/a | n/a | not planned |
| `csq` | **none** | n/a | n/a | not planned |
| `convert` | **none** | n/a | n/a | partially covered by `view -O b` and `view -O v` |
| `consensus` | **none** | n/a | n/a | not planned |
| `cnv` | **none** | n/a | n/a | not planned |

## Field-level coverage within supported subcommands

The "100 %" numbers above refer to a canonical 7-tuple
`(CHROM, POS, REF, ALT, ID, QUAL, FILTER)` plus the scalar INFO/FORMAT
fields named in the `tests/test_bcftools_*.py` suite. For exhaustive field-level
behaviour (Number=A/R/G, missing values, BND, phased GT), see
`docs/vcf_conformance_matrix.md`, which is the authoritative source.

## What rubam-bcftools does NOT do

1. **Variant calling** (`bcftools call`, `bcftools mpileup --output-VCF`).
   Out of scope. rubam is an I/O library, not a caller.
2. **Set operations across multiple VCFs** (`isec`, `merge`, `concat`).
   Roadmap v0.5.
3. **Annotation engines** (`annotate`, `csq`). Not planned.
4. **Consensus / haplotype** (`consensus`, `cnv`). Not planned.
5. **Statistics & QC** (`stats`, `roh`, `gtcheck`). Roadmap v0.5.

## Reproducing this matrix

```bash
# From the rubam repo root, with rubam installed in your venv and bcftools on PATH.
# Python subcommand parity (view/norm/concat/query/index/sort/stats):
python -m pytest tests/test_bcftools_*.py -v
# Rust-side shadow-CLI cross-check against system bcftools:
cargo test --release --no-default-features --test cargo_bcftools_cli
```

## Maintenance rule

If a row's `rubam status` changes from `none` → `partial` → `full`, a
corresponding test in the `tests/test_bcftools_*.py` suite must be added or
updated **in the same PR**. The manuscript figure that quotes any
bcftools-equivalence numbers must reference this matrix and cite the
exact commit hash.
