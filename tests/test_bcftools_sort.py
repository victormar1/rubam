"""Smoke + functional tests for rubam.tools.bcftools.sort."""
from pathlib import Path

import pytest

import rubam
from rubam.tools import bcftools


_UNSORTED_VCF = """\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=10000>
##contig=<ID=chr2,length=10000>
##contig=<ID=chr10,length=10000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr10\t100\t.\tA\tG\t30\tPASS\t.
chr1\t500\t.\tC\tT\t30\tPASS\t.
chr2\t200\t.\tG\tA\t30\tPASS\t.
chr1\t100\t.\tA\tT\t30\tPASS\t.
chr10\t50\t.\tT\tC\t30\tPASS\t.
"""


@pytest.fixture
def unsorted_vcf(tmp_path):
    p = tmp_path / "unsorted.vcf"
    p.write_text(_UNSORTED_VCF)
    return p


def test_sort_basic(unsorted_vcf, tmp_path):
    out = tmp_path / "sorted.vcf"
    n = bcftools.sort(str(unsorted_vcf), str(out))
    assert n == 5
    with rubam.VariantFile(str(out), "r") as f:
        rows = [(r.reference_name, r.position) for r in f]
    # chr1 before chr2 before chr10 (header-declared order, not alphabetical)
    assert rows == [
        ("chr1", 100),
        ("chr1", 500),
        ("chr2", 200),
        ("chr10", 50),
        ("chr10", 100),
    ]


def test_sort_already_sorted_idempotent(tmp_path):
    src = tmp_path / "src.vcf"
    src.write_text("""\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30\tPASS\t.
chr1\t200\t.\tC\tT\t30\tPASS\t.
chr1\t300\t.\tG\tA\t30\tPASS\t.
""")
    out = tmp_path / "out.vcf"
    n = bcftools.sort(str(src), str(out))
    assert n == 3
    with rubam.VariantFile(str(out), "r") as f:
        positions = [r.position for r in f]
    assert positions == [100, 200, 300]


def test_sort_output_bgzf(unsorted_vcf, tmp_path):
    out = tmp_path / "sorted.vcf.gz"
    n = bcftools.sort(str(unsorted_vcf), str(out), output_type="z")
    assert n == 5
    # First 2 bytes are gzip magic
    assert out.read_bytes()[:2] == b"\x1f\x8b"


def test_sort_output_bcf(unsorted_vcf, tmp_path):
    out = tmp_path / "sorted.bcf"
    n = bcftools.sort(str(unsorted_vcf), str(out), output_type="b")
    assert n == 5
    # Round-trip: read it back and check order
    with rubam.VariantFile(str(out), "r") as f:
        rows = [(r.reference_name, r.position) for r in f]
    assert len(rows) == 5
    assert rows[0] == ("chr1", 100)


def test_sort_unknown_format_raises(unsorted_vcf, tmp_path):
    with pytest.raises(ValueError):
        bcftools.sort(str(unsorted_vcf), str(tmp_path / "x.vcf"), output_type="q")


def test_sort_empty_vcf(tmp_path):
    src = tmp_path / "empty.vcf"
    src.write_text("""\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
""")
    out = tmp_path / "out.vcf"
    n = bcftools.sort(str(src), str(out))
    assert n == 0
