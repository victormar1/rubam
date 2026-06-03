"""Smoke + functional tests for rubam.tools.bcftools.view."""
import shutil
import subprocess
from pathlib import Path

import pytest

import rubam
from rubam.tools import bcftools


_VCF_TEXT = """##fileformat=VCFv4.3
##contig=<ID=chr1,length=10000>
##contig=<ID=chr2,length=10000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\tNA12891\tNA12892
chr1\t100\t.\tA\tG\t30\tPASS\tDP=20\tGT\t0/1\t1/1\t0/0
chr1\t500\t.\tC\tT\t40\tPASS\tDP=25\tGT\t0/0\t0/1\t1/1
chr2\t200\t.\tG\tA\t50\tPASS\tDP=30\tGT\t0/1\t1/1\t0/1
"""


@pytest.fixture
def src_vcf(tmp_path):
    p = tmp_path / "src.vcf"
    p.write_text(_VCF_TEXT)
    return p


def test_view_count_only_no_filter(src_vcf, tmp_path):
    out = tmp_path / "out.vcf"
    n = bcftools.view(str(src_vcf), output=str(out))
    assert n == 3
    text = out.read_text()
    assert "chr1\t100" in text
    assert "chr2\t200" in text


def test_view_region_chr1(src_vcf, tmp_path):
    out = tmp_path / "out.vcf"
    n = bcftools.view(str(src_vcf), region="chr1:1-1000", output=str(out))
    assert n == 2  # chr1 100 and chr1 500
    text = out.read_text()
    assert "chr1\t100" in text
    assert "chr1\t500" in text
    assert "chr2" not in text or "##contig=<ID=chr2" in text  # contig in header is OK


def test_view_region_chr2_only(src_vcf, tmp_path):
    out = tmp_path / "out.vcf"
    n = bcftools.view(str(src_vcf), region="chr2", output=str(out))
    assert n == 1
    text = out.read_text()
    assert "chr2\t200" in text


def test_view_subset_samples(src_vcf, tmp_path):
    out = tmp_path / "out.vcf"
    n = bcftools.view(
        str(src_vcf), samples=["NA12878", "NA12892"], output=str(out)
    )
    assert n == 3
    # Check the output header has only the kept samples in the column order
    with rubam.VariantFile(str(out), "r") as f:
        assert f.header.samples == ("NA12878", "NA12892")


def test_view_header_only(src_vcf, tmp_path):
    out = tmp_path / "out.vcf"
    n = bcftools.view(str(src_vcf), header_only=True, output=str(out))
    assert n == 0
    text = out.read_text()
    assert "##fileformat=" in text
    assert "chr1\t100" not in text


def test_view_no_header(src_vcf, tmp_path):
    out = tmp_path / "out.vcf"
    n = bcftools.view(str(src_vcf), no_header=True, output=str(out))
    assert n == 3
    text = out.read_text()
    assert "##fileformat=" not in text
    assert "chr1\t100" in text


def test_view_unknown_format_raises(src_vcf, tmp_path):
    with pytest.raises(ValueError):
        bcftools.view(str(src_vcf), output_type="q", output=str(tmp_path / "x.vcf"))


def test_view_output_bgzf(src_vcf, tmp_path):
    out = tmp_path / "out.vcf.gz"
    n = bcftools.view(str(src_vcf), output_type="z", output=str(out))
    assert n == 3
    # File should be BGZF — first 2 bytes are 0x1f 0x8b (gzip magic)
    data = out.read_bytes()
    assert data[:2] == b"\x1f\x8b"
    # Read it back to confirm content survives
    with rubam.VariantFile(str(out), "r") as f:
        assert sum(1 for _ in f) == 3
