"""rubam command-line interface.

Sub-commands:
    rubam depth        Per-base depth over a region.
    rubam depth-bed    Per-base depth for every region in a BED file.
    rubam pileup       A/C/G/T/N counts per position over a region.
    rubam count        Count reads matching SAM-flag and MAPQ filters.
    rubam flagstat     samtools flagstat replacement (JSON or text).
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Iterable

from rubam import (
    __version__,
    count_reads,
    flag_stats,
    get_depths,
    get_depths_regions,
    pileup_bases,
)


# ---------- shared option helpers ---------- #

def _add_filter_options(p: argparse.ArgumentParser) -> None:
    p.add_argument("-Q", "--min-mapq", type=int, default=0, dest="min_mapq",
                   help="Minimum mapping quality (default: %(default)s)")
    p.add_argument("-q", "--min-bq", type=int, default=13, dest="min_bq",
                   help="Minimum base quality (default: %(default)s)")
    p.add_argument("-d", "--max-depth", type=int, default=8000, dest="max_depth",
                   help="Cap on reported depth per position (default: %(default)s)")
    p.add_argument("-n", "--num-threads", type=int, default=12, dest="num_threads",
                   help="Worker threads (default: %(default)s)")
    p.add_argument("-t", "--step", type=int, default=1,
                   help="Sampling step (default: %(default)s)")


def _read_bed(path: str) -> list[tuple[str, int, int]]:
    out: list[tuple[str, int, int]] = []
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                raise ValueError(f"BED line has fewer than 3 columns: {raw!r}")
            chrom = parts[0]
            # BED is 0-based half-open; rubam APIs are 1-based inclusive.
            start = int(parts[1]) + 1
            end = int(parts[2])
            if end < start:
                continue
            out.append((chrom, start, end))
    return out


# ---------- sub-command implementations ---------- #

def _add_depth_subparser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("depth",
                       help="Per-base sequencing depth over a region.")
    p.add_argument("bam")
    p.add_argument("chromosome")
    p.add_argument("start", type=int)
    p.add_argument("end", type=int)
    _add_filter_options(p)
    p.add_argument("-j", "--json", action="store_true",
                   help="Emit JSON instead of TSV")
    p.set_defaults(func=_cmd_depth)


def _cmd_depth(args: argparse.Namespace) -> int:
    positions, depths = get_depths(
        args.bam, args.chromosome, args.start, args.end,
        args.step, args.min_mapq, args.min_bq, args.max_depth, args.num_threads,
    )
    if args.json:
        json.dump(dict(zip(positions, depths)), sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        for pos, depth in zip(positions, depths):
            sys.stdout.write(f"{args.chromosome}\t{pos}\t{depth}\n")
    return 0


def _add_depth_bed_subparser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("depth-bed",
                       help="Per-base depth for every region in a BED file.")
    p.add_argument("bam")
    p.add_argument("bed")
    _add_filter_options(p)
    p.set_defaults(func=_cmd_depth_bed)


def _cmd_depth_bed(args: argparse.Namespace) -> int:
    regions = _read_bed(args.bed)
    results = get_depths_regions(
        args.bam, regions,
        args.step, args.min_mapq, args.min_bq, args.max_depth, args.num_threads,
    )
    out = sys.stdout
    for (chrom, _, _), (positions, depths) in zip(regions, results):
        for pos, depth in zip(positions, depths):
            out.write(f"{chrom}\t{pos}\t{depth}\n")
    return 0


def _add_pileup_subparser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("pileup",
                       help="A/C/G/T/N counts per position over a region.")
    p.add_argument("bam")
    p.add_argument("chromosome")
    p.add_argument("start", type=int)
    p.add_argument("end", type=int)
    _add_filter_options(p)
    p.set_defaults(func=_cmd_pileup)


def _cmd_pileup(args: argparse.Namespace) -> int:
    positions, a, c, g, t, n, depth = pileup_bases(
        args.bam, args.chromosome, args.start, args.end,
        args.step, args.min_mapq, args.min_bq, args.max_depth, args.num_threads,
    )
    out = sys.stdout
    out.write("chrom\tpos\tdepth\tA\tC\tG\tT\tN\n")
    for i, pos in enumerate(positions):
        out.write(
            f"{args.chromosome}\t{pos}\t{depth[i]}\t"
            f"{a[i]}\t{c[i]}\t{g[i]}\t{t[i]}\t{n[i]}\n"
        )
    return 0


def _add_count_subparser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("count",
                       help="Count reads matching SAM-flag and MAPQ filters.")
    p.add_argument("bam")
    p.add_argument("chromosome")
    p.add_argument("start", type=int)
    p.add_argument("end", type=int)
    p.add_argument("-Q", "--min-mapq", type=int, default=0, dest="min_mapq")
    p.add_argument("-f", "--flag-required", type=lambda s: int(s, 0),
                   default=0, dest="flag_required",
                   help="Required flags (e.g. 0x2, default 0)")
    p.add_argument("-F", "--flag-filtered", type=lambda s: int(s, 0),
                   default=0x704, dest="flag_filtered",
                   help="Excluded flags (default 0x704: UNMAP|SECONDARY|QCFAIL|DUP)")
    p.set_defaults(func=_cmd_count)


def _cmd_count(args: argparse.Namespace) -> int:
    n = count_reads(
        args.bam, args.chromosome, args.start, args.end,
        args.min_mapq, args.flag_required, args.flag_filtered,
    )
    sys.stdout.write(f"{n}\n")
    return 0


def _add_flagstat_subparser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("flagstat",
                       help="samtools flagstat replacement.")
    p.add_argument("bam")
    p.add_argument("-j", "--json", action="store_true",
                   help="Emit JSON instead of human text")
    p.set_defaults(func=_cmd_flagstat)


_FLAGSTAT_ORDER = [
    ("total", "in total (QC-passed reads)"),
    ("qcfail", "QC-failed reads"),
    ("primary", "primary"),
    ("secondary", "secondary"),
    ("supplementary", "supplementary"),
    ("duplicates", "duplicates"),
    ("primary_duplicates", "primary duplicates"),
    ("mapped", "mapped"),
    ("primary_mapped", "primary mapped"),
    ("paired", "paired in sequencing"),
    ("read_1", "read1"),
    ("read_2", "read2"),
    ("properly_paired", "properly paired"),
    ("with_itself_and_mate_mapped", "with itself and mate mapped"),
    ("singletons", "singletons"),
    ("mate_mapped_to_different_chr", "with mate mapped to a different chr"),
    ("mate_mapped_to_different_chr_mapq_5",
     "with mate mapped to a different chr (mapQ>=5)"),
]


def _cmd_flagstat(args: argparse.Namespace) -> int:
    stats = flag_stats(args.bam)
    if args.json:
        json.dump(stats, sys.stdout, indent=2)
        sys.stdout.write("\n")
        return 0
    for key, label in _FLAGSTAT_ORDER:
        sys.stdout.write(f"{stats.get(key, 0)}\t{label}\n")
    return 0


# ---------- entry point ---------- #

def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="rubam",
        description="Pure-Rust BAM/CRAM coverage, pileup and stats. "
                    "Native Windows / Linux / macOS, multi-threaded.",
    )
    parser.add_argument("--version", action="version",
                        version=f"rubam {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)
    _add_depth_subparser(sub)
    _add_depth_bed_subparser(sub)
    _add_pileup_subparser(sub)
    _add_count_subparser(sub)
    _add_flagstat_subparser(sub)

    args = parser.parse_args(argv if argv is None else list(argv))
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
