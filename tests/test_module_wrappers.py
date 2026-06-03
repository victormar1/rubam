"""Tests for the pysam-compatible module-level wrappers exposed on the
top-level `rubam` namespace: rubam.flagstat, rubam.view, rubam.idxstats,
rubam.faidx, rubam.depth, rubam.merge.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

import rubam


EXAMPLE_BAM = str(Path(__file__).parent / "example.bam")

_TINY_FA = ">chr1\nACGTACGTACGTACGTACGT\n>chr2\nNNNNNNNNNN\n"


@pytest.fixture
def tiny_fa(tmp_path: Path) -> str:
    fa = tmp_path / "tiny.fa"
    fa.write_text(_TINY_FA)
    return str(fa)


def test_module_flagstat_returns_multiline_string():
    out = rubam.flagstat(EXAMPLE_BAM)
    assert isinstance(out, str)
    lines = out.strip().splitlines()
    # samtools flagstat is 13-17 lines depending on version; we accept >= 10.
    assert len(lines) >= 10
    # Every line is "<n> + 0 <label>"
    assert "in total" in out
    assert "mapped" in out


def test_module_idxstats_returns_list_of_dicts():
    out = rubam.idxstats(EXAMPLE_BAM)
    assert isinstance(out, list)
    assert len(out) >= 1
    row = out[0]
    assert "contig" in row
    assert "length" in row
    assert "mapped" in row
    assert "unmapped" in row


def test_module_view_count_only_returns_int():
    n = rubam.view(EXAMPLE_BAM, count_only=True)
    assert isinstance(n, int)
    assert n > 0


def test_module_view_returns_record_strings():
    """Default view (no output= and not count_only) returns list of str."""
    out = rubam.view(EXAMPLE_BAM)
    assert isinstance(out, list)
    assert len(out) > 0
    for entry in out[:5]:
        assert isinstance(entry, str)


def test_module_view_with_region():
    """rubam.view(bam, region) should restrict to the region."""
    n_all = rubam.view(EXAMPLE_BAM, count_only=True)
    n_region = rubam.view(EXAMPLE_BAM, region="chr1:1-1000000", count_only=True)
    assert n_region <= n_all


def test_module_faidx_no_region_builds_index(tmp_path, tiny_fa):
    # The fixture already has no .fai; faidx() builds it.
    fai = tiny_fa + ".fai"
    if os.path.exists(fai):
        os.remove(fai)
    result = rubam.faidx(tiny_fa)
    assert result is None
    assert os.path.exists(fai)


def test_module_faidx_subseq(tiny_fa):
    rubam.faidx(tiny_fa)  # build .fai
    out = rubam.faidx(tiny_fa, "chr1:1-10")
    # pysam.faidx returns FASTA-formatted output: ">chr1:1-10\nACGTACGTAC\n"
    assert isinstance(out, str)
    assert ">chr1:1-10" in out
    assert "ACGTACGTAC" in out


def test_module_faidx_multiple_regions(tiny_fa):
    rubam.faidx(tiny_fa)
    out = rubam.faidx(tiny_fa, "chr1:1-4", "chr2:1-5")
    assert "ACGT" in out
    assert "NNNNN" in out


def test_module_depth_returns_tsv_string():
    out = rubam.depth(EXAMPLE_BAM, region="chr1:1000000-1000050")
    assert isinstance(out, str)
    if not out.strip():
        # zero coverage in that region — still a valid (empty) TSV
        return
    for line in out.strip().splitlines():
        cols = line.split("\t")
        assert len(cols) == 3
        # Second and third cols parse as ints
        int(cols[1])
        int(cols[2])


def test_module_depth_explicit_coords():
    out = rubam.depth(EXAMPLE_BAM, chromosome="chr1", start=1000000, end=1000010)
    assert isinstance(out, str)


def test_module_depth_invalid_args_raise():
    with pytest.raises(ValueError):
        rubam.depth(EXAMPLE_BAM)  # no region or coords


def test_module_merge_concatenates_bams(tmp_path):
    # Merge example.bam with itself — exercises the call shape.
    out = tmp_path / "merged.bam"
    rubam.merge(str(out), EXAMPLE_BAM, EXAMPLE_BAM)
    assert out.exists()
    assert out.stat().st_size > 0


def test_module_merge_requires_inputs(tmp_path):
    out = tmp_path / "merged.bam"
    with pytest.raises(ValueError):
        rubam.merge(str(out))
