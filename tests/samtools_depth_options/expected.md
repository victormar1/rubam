# `samtools depth` option compatibility — empirical matrix

**Generated:** 2026-05-13 against rubam release binaries (Windows
`target/release/rubam-samtools.exe`, `target/release/rubam-depth.exe`)
on WSL Ubuntu 22.04 with system samtools 1.13 (`apt`).

> **Headline finding:** `rubam-samtools depth` is **not implemented**
> as a subcommand. Every option below returns
> `rubam-samtools: unknown subcommand "depth"` from the shadow CLI.
> The depth functionality is shipped as a separate binary,
> `rubam-depth`, with its own argv surface. The "samtools-compatible"
> manuscript claim must therefore be qualified to the `rubam-depth`
> binary, not the `rubam-samtools depth` shape.

Fixture: `data/synthetic/synth_chr20_1Mb_30x.bam` (9.8 MB, single
chromosome `chr20`, 1 Mb, 30x synthetic coverage, MAPQ 60, BQ 30).
Region for every test: `chr20:1-1000`.

Reproduction:

```bash
bash tests/samtools_depth_options/run_matrix.sh
```

## Status legend

| status        | meaning                                                             |
|---------------|---------------------------------------------------------------------|
| **yes**       | `rubam-samtools depth <args>` accepts the option and matches samtools |
| **partial**   | accepted but output diverges                                         |
| **no**        | not accepted by `rubam-samtools depth`                              |
| *fallback*    | the standalone `rubam-depth` binary may cover the same intent       |

Output match column applies to the **rubam-samtools** invocation
unless explicitly marked "(via rubam-depth)".

## Matrix

| option       | supported | output match           | notes                                                                                                                                  |
|--------------|-----------|------------------------|----------------------------------------------------------------------------------------------------------------------------------------|
| *(no flag)*  | no        | n/a                    | `unknown subcommand "depth"`. **Default behaviour diverges**: samtools skips zero-depth positions, `rubam-depth` always emits them.   |
| `-a`         | no        | n/a                    | `unknown subcommand "depth"`. *Fallback:* `rubam-depth` output is **exact** on `chr20:1-1000` (`-a` matches rubam-depth's default).    |
| `-aa`        | no        | n/a                    | `unknown subcommand "depth"`. *Fallback:* `rubam-depth` is exact within the requested region; unused-ref behaviour untested (fixture has one chromosome). |
| `-q 13`      | no        | n/a                    | `unknown subcommand "depth"`. *Fallback:* `rubam-depth -q 13` is **exact** vs `samtools depth -a -q 13`. Note `rubam-depth`'s default is `-q 13`, samtools's is `-q 0`. |
| `-Q 10`      | no        | n/a                    | `unknown subcommand "depth"`. *Fallback:* `rubam-depth -Q 10` is **exact** vs `samtools depth -a -Q 10`.                                |
| `-r REGION`  | no        | n/a                    | `unknown subcommand "depth"`. *Fallback:* `rubam-depth` takes the region as 3 positional args (`<chrom> <start> <end>`) — same semantics, different shape. |
| `-b BED`     | no        | n/a                    | `unknown subcommand "depth"`. **No `rubam-depth` equivalent**: only single-region invocation is supported. Manuscript must drop the BED claim. |
| `-G 0x4`     | no        | n/a                    | `unknown subcommand "depth"`. **No `rubam-depth` equivalent**: no flag-filter knob; UNMAP/SECONDARY/etc. filtering not exposed.        |
| `-d 100`     | no        | n/a                    | `unknown subcommand "depth"`. *Fallback:* `rubam-depth -d 100` is **exact** vs `samtools depth -a -d 100`.                              |
| `-H`         | no        | n/a                    | `unknown subcommand "depth"`. **No `rubam-depth` equivalent**: rubam never emits a header row.                                         |
| `-o FILE`    | no        | n/a                    | `unknown subcommand "depth"`. **No `rubam-depth` equivalent**: output is stdout only (the caller pipes/redirects).                     |

## Summary

| count | category                                                              |
|-------|-----------------------------------------------------------------------|
| 0     | options supported by `rubam-samtools depth`                           |
| 0     | options partially supported by `rubam-samtools depth`                 |
| 11    | options not implemented (incl. default behaviour)                     |
| 6     | options covered with **exact** output by the `rubam-depth` fallback   |
| 4     | options with **no** `rubam-depth` equivalent (`-b`, `-G`, `-H`, `-o`) |
| 1     | options where `rubam-depth` **diverges** from samtools default (no flag → emits zero-depth rows samtools omits) |

## Reproduction artefacts

- Per-option stdout/stderr captures: `tests/samtools_depth_options/results/`
- Machine-readable summary TSV: `tests/samtools_depth_options/results/summary.tsv`
- BED file used for `-b`: `tests/samtools_depth_options/results/region.bed`

## Recommended manuscript phrasing

Replace any claim of *"drop-in `samtools depth` replacement"* with:

> rubam ships a standalone `rubam-depth` CLI that reproduces
> `samtools depth -a` output byte-exactly on the canonical option
> subset `{-a, -aa (within-region), -q, -Q, -r, -d}` for indexed BAMs.
> The `-b`, `-G`, `-H`, `-o`, and multi-BAM input options of
> `samtools depth` are **not** implemented; the depth subcommand of
> the wider `rubam-samtools` shadow CLI is **not** wired in this
> release.
