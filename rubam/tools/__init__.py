"""Pure-Rust ports of samtools / bcftools subcommands.

This package wraps the Rust implementations in src/tools/. Each function takes
the same arguments as the corresponding samtools subcommand whenever possible.
Implementations are filled in across tasks B2..B8.
"""

from __future__ import annotations

from rubam._rubam import calmd as _calmd
from rubam._rubam import faidx as _faidx
from rubam._rubam import flag_stats as _flag_stats_v01
from rubam._rubam import idxstats as _idxstats
from rubam._rubam import index as _index
from rubam._rubam import merge as _merge
from rubam._rubam import sort as _sort
from rubam._rubam import view as _view

from rubam.tools import bcftools  # noqa: F401  re-export bcftools sub-namespace


def flagstat(input):
    """samtools flagstat equivalent. Returns a dict.

    Thin forwarder over rubam.flag_stats; the two return identical dicts.
    """
    return _flag_stats_v01(input)


def sort(input: str, output: str, *, threads: int = 1) -> None:
    """Coordinate-sort a BAM. Equivalent to `samtools sort INPUT -o OUTPUT`."""
    _sort(input, output, threads=threads)


def index(input: str, *, csi: bool = False) -> None:
    """Build a .bai index for a coordinate-sorted BAM. Equivalent to `samtools index`."""
    _index(input, csi=csi)


def view(input: str, *, region=None, output=None,
         min_mapq=0, flag_required=0, flag_filtered=0, count_only=False):
    """samtools view equivalent.

    If count_only=True, returns the matching record count and ignores output.
    Otherwise, if output is set, writes a filtered BAM and still returns the count.
    """
    return _view(
        input, region=region, output=output,
        min_mapq=min_mapq, flag_required=flag_required, flag_filtered=flag_filtered,
        count_only=count_only,
    )


def merge(inputs, output, *, force=True):
    """Merge multiple coordinate-sorted BAMs into one. Equivalent to `samtools merge`."""
    _merge(list(inputs), output, force=force)


def idxstats(input: str):
    """samtools idxstats equivalent. Returns a list of dicts with
    keys 'contig', 'length', 'mapped', 'unmapped'."""
    return _idxstats(input)


def faidx(input: str, *, region: str | None = None):
    """samtools faidx equivalent. Writes .fai if missing.

    If region (chr:start-end, 1-based inclusive) is given, returns the
    subsequence as a string. Otherwise returns None.
    """
    return _faidx(input, region=region)


def calmd(input: str, reference: str, *, output: str):
    """samtools calmd port. v0.2 emits NM only; MD lands in v0.2.x."""
    _calmd(input, reference, output=output)


__all__ = ["calmd", "faidx", "flagstat", "idxstats", "index", "merge", "sort", "view"]
