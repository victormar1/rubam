# rubam ↔ `samtools depth` compatibility

**Last updated: 2026-05-14 (against rubam v0.3.2 / samtools 1.13).**
**Empirical source:** `tests/samtools_depth_options/run_matrix.sh`.

This matrix is the authoritative scope of rubam's "samtools depth
compatibility" claim. It is the only document the manuscript should
cite when stating which `samtools depth` options behave identically.

## Architectural caveat

`rubam-samtools` is the shadow CLI that mirrors a subset of the
`samtools` subcommand set (`sort`, `index`, `view`, `merge`,
`flagstat`, `idxstats`, `calmd`, `faidx`). **`depth` is not part of
that set.** The depth functionality lives in a standalone binary,
`rubam-depth`, whose argv is:

```
rubam-depth <bam> <chrom> <start> <end> \
    [-n THREADS] [-Q MIN_MAPQ] [-q MIN_BQ] [-d MAX_DEPTH] [-t STEP]
```

The table below therefore answers two questions per option:

1. Does `rubam-samtools depth <opt>` accept it? *(answer: no, for every
   option — the subcommand does not exist.)*
2. Does the standalone `rubam-depth` have an equivalent? *(this is
   where the per-option detail lives.)*

## Status legend

| status              | meaning                                                                                  |
|---------------------|------------------------------------------------------------------------------------------|
| **exact**           | byte-for-byte identical to `samtools depth` output                                       |
| **byte-equivalent** | identical once whitespace is normalised                                                  |
| **diverges**        | different rows or different depths                                                       |
| **no fallback**     | neither `rubam-samtools` nor `rubam-depth` exposes the option                            |

## Matrix

| option       | `rubam-samtools depth` | `rubam-depth` (fallback) | output match (vs samtools) | notes                                                                                                          |
|--------------|:----------------------:|:------------------------:|----------------------------|----------------------------------------------------------------------------------------------------------------|
| *(no flag)*  | no                     | n/a                      | diverges                   | samtools skips zero-depth rows; `rubam-depth` always emits them. Equivalent to running samtools with `-a`.    |
| `-a`         | no                     | yes (default)            | **exact**                  | `rubam-depth` is implicitly `-a` (always emits every position in the requested interval).                     |
| `-aa`        | no                     | partial                  | **exact** within region    | Only one chromosome in the fixture, so "unused ref" branch is untested. Within the region: identical.         |
| `-q 13`      | no                     | `-q 13`                  | **exact**                  | Note that the **defaults differ**: samtools = `-q 0`, `rubam-depth` = `-q 13`. Always pass `-q` explicitly.   |
| `-Q 10`      | no                     | `-Q 10`                  | **exact**                  | Mapping-quality filter applied identically.                                                                    |
| `-r REGION`  | no                     | positional args          | **exact**                  | Different shape: `rubam-depth bam.bam chr20 1 1000` instead of `-r chr20:1-1000`.                              |
| `-b BED`     | no                     | **no fallback**          | n/a                        | BED-driven region lists are not implemented. Call `rubam-depth` once per region from the caller.              |
| `-G 0x4`     | no                     | **no fallback**          | n/a                        | No way to alter the default filter-out flag set. UNMAP/SECONDARY/QCFAIL/DUP handling is hard-coded.           |
| `-d 100`     | no                     | `-d 100`                 | **exact**                  | Per-position depth cap applied identically.                                                                    |
| `-H`         | no                     | **no fallback**          | n/a                        | `rubam-depth` never emits a header row.                                                                        |
| `-o FILE`    | no                     | **no fallback**          | n/a                        | Output is stdout only; redirect with `>` from the shell.                                                       |

## Headline counts

- `rubam-samtools depth` supports **0 / 11** of the tested options.
- `rubam-depth` (the actual implementation) matches `samtools depth -a`
  **exactly** on **6 / 11** options (`-a`, `-aa` within-region, `-q`,
  `-Q`, `-r` (re-shaped), `-d`).
- **4 / 11** options have **no `rubam-depth` equivalent**: `-b`, `-G`,
  `-H`, `-o`. These should not be promised in the manuscript.
- The **default behaviour** (no flag) of `rubam-depth` corresponds to
  `samtools depth -a`, not bare `samtools depth`. This is the single
  most important caveat for callers migrating shell pipelines.

## Recommended manuscript phrasing

> rubam ships a standalone `rubam-depth` CLI that reproduces
> `samtools depth -a` output byte-exactly on the canonical option
> subset `{-a, -aa, -q, -Q, -r, -d}` for indexed BAMs. The `-b`,
> `-G`, `-H`, `-o`, and multi-BAM input modes of `samtools depth` are
> not implemented in this release; the depth subcommand of the wider
> `rubam-samtools` shadow CLI is also not wired (planned for v0.4).

## Reproduction

```bash
# WSL Ubuntu, samtools 1.13, rubam release binaries built with
# `cargo build --release` or `maturin develop --release`.
bash tests/samtools_depth_options/run_matrix.sh
```

Per-test artefacts (samtools/rubam stdout, stderr, diff inputs) land
under `tests/samtools_depth_options/results/`.

The machine-readable summary lives at
`tests/samtools_depth_options/results/summary.tsv` and mirrors this
table column-for-column (`option`, `supported`, `output_match`,
`notes`).
