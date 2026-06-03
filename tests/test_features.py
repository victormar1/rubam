"""Smoke tests for the v0.1 feature surface beyond `get_depths`.

The fixtures rely on the same `tests/example.bam` used by `test_core.py`. The
expected golden values are taken from the existing depth fixture so the new
APIs are exercised against known-good ground truth.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from rubam import (
    count_reads,
    flag_stats,
    get_depths,
    get_depths_regions,
    pileup_bases,
)


EXAMPLE_BAM = str(Path(__file__).parent / "example.bam")
CHROM = "chr1"
START = 1_000_000
END = 1_000_020

# Same golden values as `test_core.py::test_get_depths`.
GOLDEN_DEPTHS = [
    51, 52, 44, 52, 53, 47, 51, 52, 49, 50,
    49, 50, 50, 49, 50, 50, 46, 50, 48, 50, 44,
]


def test_count_reads_matches_depth_floor():
    """`count_reads` over the region must be >= max depth in that region.

    Every read counted at one position must have been observed by the count;
    therefore `count_reads >= max(depths)`. We don't require equality because a
    read can span multiple positions but is counted once by `count_reads`.
    """
    n = count_reads(EXAMPLE_BAM, CHROM, START, END)
    assert n >= max(GOLDEN_DEPTHS), n


def test_count_reads_default_excludes_dup_secondary():
    """Default flag mask = UNMAP|SECONDARY|QCFAIL|DUP. Asking for the
    *complement* (i.e. requiring those flags to be set) should yield 0 here."""
    n = count_reads(EXAMPLE_BAM, CHROM, START, END,
                    flag_required=0x100, flag_filtered=0)
    assert n == 0, n


def test_pileup_bases_sums_to_depth():
    positions, a, c, g, t, n, depth = pileup_bases(
        EXAMPLE_BAM, CHROM, START, END,
    )
    assert positions == list(range(START, END + 1))
    assert depth == GOLDEN_DEPTHS, depth
    for i in range(len(positions)):
        assert a[i] + c[i] + g[i] + t[i] + n[i] == depth[i], (
            i, a[i], c[i], g[i], t[i], n[i], depth[i],
        )


def test_pileup_bases_all_threads_match():
    """The chunked parallel implementation must be deterministic w.r.t. threads."""
    ref = pileup_bases(EXAMPLE_BAM, CHROM, START, END, num_threads=1)
    for nt in (2, 4, 8, 16):
        out = pileup_bases(EXAMPLE_BAM, CHROM, START, END, num_threads=nt)
        assert out == ref, nt


def test_get_depths_regions_matches_individual_calls():
    regions = [
        (CHROM, START, START + 5),
        (CHROM, START + 10, END),
    ]
    batch = get_depths_regions(EXAMPLE_BAM, regions)
    individual = [
        get_depths(EXAMPLE_BAM, c, s, e) for (c, s, e) in regions
    ]
    assert batch == individual


def test_get_depths_regions_handles_empty_list():
    assert get_depths_regions(EXAMPLE_BAM, []) == []


def test_flag_stats_returns_dict_with_expected_keys():
    stats = flag_stats(EXAMPLE_BAM)
    assert isinstance(stats, dict)
    expected_keys = {
        "total", "qcfail", "primary", "secondary", "supplementary",
        "duplicates", "primary_duplicates", "mapped", "primary_mapped",
        "paired", "read_1", "read_2", "properly_paired",
        "with_itself_and_mate_mapped", "singletons",
        "mate_mapped_to_different_chr", "mate_mapped_to_different_chr_mapq_5",
    }
    assert expected_keys.issubset(stats.keys())
    # Sanity invariants.
    assert stats["total"] > 0
    assert stats["primary"] <= stats["total"]
    assert stats["primary_mapped"] <= stats["primary"]
    assert stats["mapped"] >= stats["primary_mapped"]
    assert stats["paired"] >= stats["read_1"]
    assert stats["paired"] >= stats["read_2"]


def test_flag_stats_invalid_path_raises():
    with pytest.raises((IOError, OSError)):
        flag_stats("does_not_exist_zzz.bam")
