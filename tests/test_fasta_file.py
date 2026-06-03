"""Tests for `rubam.FastaFile` — pysam.FastaFile drop-in."""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pytest

import rubam


# A tiny FASTA fixture (built inline so the test runs without external data).
_TINY_FA = ">chr1\nACGTACGTACGTACGTACGT\n>chr2\nNNNNNNNNNN\n"


@pytest.fixture
def tiny_fa(tmp_path: Path) -> str:
    fa = tmp_path / "tiny.fa"
    fa.write_text(_TINY_FA)
    return str(fa)


def test_open_and_close(tiny_fa: str):
    fa = rubam.FastaFile(tiny_fa)
    assert fa.is_open
    assert fa.references == ("chr1", "chr2")
    assert fa.lengths == (20, 10)
    assert fa.nreferences == 2
    fa.close()
    assert not fa.is_open


def test_context_manager(tiny_fa: str):
    with rubam.FastaFile(tiny_fa) as fa:
        assert fa.is_open
        assert fa.fetch("chr1") == "ACGTACGTACGTACGTACGT"


def test_fetch_full_contig(tiny_fa: str):
    with rubam.FastaFile(tiny_fa) as fa:
        assert fa.fetch("chr1") == "ACGTACGTACGTACGTACGT"
        assert fa.fetch("chr2") == "NNNNNNNNNN"


def test_fetch_slice_zero_based_halfopen(tiny_fa: str):
    """pysam convention: fetch('chr1', 0, 4) returns the first 4 bases."""
    with rubam.FastaFile(tiny_fa) as fa:
        assert fa.fetch("chr1", 0, 4) == "ACGT"
        assert fa.fetch("chr1", 4, 8) == "ACGT"
        assert fa.fetch("chr1", 0, 1) == "A"


def test_fetch_empty_range(tiny_fa: str):
    with rubam.FastaFile(tiny_fa) as fa:
        assert fa.fetch("chr1", 5, 5) == ""


def test_fetch_region_string(tiny_fa: str):
    """samtools 1-based inclusive region syntax via region= kwarg."""
    with rubam.FastaFile(tiny_fa) as fa:
        # chr1:1-4 is 1-based inclusive = bases 1-4 = "ACGT".
        assert fa.fetch(region="chr1:1-4") == "ACGT"


def test_fetch_unknown_contig_raises(tiny_fa: str):
    with rubam.FastaFile(tiny_fa) as fa:
        with pytest.raises(KeyError):
            fa.get_reference_length("chrUNK")
        with pytest.raises((KeyError, OSError, IOError)):
            fa.fetch("chrUNK", 0, 5)


def test_fetch_out_of_range_raises(tiny_fa: str):
    with rubam.FastaFile(tiny_fa) as fa:
        with pytest.raises(ValueError):
            fa.fetch("chr1", 0, 100)


def test_fetch_negative_raises(tiny_fa: str):
    with rubam.FastaFile(tiny_fa) as fa:
        with pytest.raises(ValueError):
            fa.fetch("chr1", -1, 5)


def test_fetch_after_close_raises(tiny_fa: str):
    fa = rubam.FastaFile(tiny_fa)
    fa.close()
    with pytest.raises((OSError, IOError)):
        fa.fetch("chr1", 0, 5)


def test_get_reference_length(tiny_fa: str):
    with rubam.FastaFile(tiny_fa) as fa:
        assert fa.get_reference_length("chr1") == 20
        assert fa.get_reference_length("chr2") == 10


def test_pathlib_path_accepted(tiny_fa: str):
    p = Path(tiny_fa)
    with rubam.FastaFile(p) as fa:
        assert fa.fetch("chr1", 0, 4) == "ACGT"


def test_index_auto_built(tmp_path: Path):
    """If the .fai is missing, FastaFile builds it on open."""
    fa_path = tmp_path / "fresh.fa"
    fa_path.write_text(_TINY_FA)
    # No .fai yet.
    assert not (tmp_path / "fresh.fa.fai").exists()
    fa = rubam.FastaFile(str(fa_path))
    assert (tmp_path / "fresh.fa.fai").exists()
    fa.close()
