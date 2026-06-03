"""Smoke + functional tests for rubam.tools.bcftools.concat."""
import pytest
import rubam
from rubam.tools import bcftools


_HEADER = """##fileformat=VCFv4.3
##contig=<ID=chr1,length=10000>
##contig=<ID=chr2,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878
"""

_PART_A = _HEADER + "chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0/1\n"
_PART_B = _HEADER + "chr1\t500\t.\tC\tT\t40\tPASS\t.\tGT\t1/1\nchr2\t200\t.\tG\tA\t50\tPASS\t.\tGT\t0/0\n"


def _write(p, txt):
    p.write_text(txt)
    return p


def test_concat_two_files(tmp_path):
    a = _write(tmp_path / "a.vcf", _PART_A)
    b = _write(tmp_path / "b.vcf", _PART_B)
    out = tmp_path / "out.vcf"
    n = bcftools.concat([str(a), str(b)], str(out))
    assert n == 3
    with rubam.VariantFile(str(out), "r") as f:
        positions = [(r.reference_name, r.position) for r in f]
    assert positions == [("chr1", 100), ("chr1", 500), ("chr2", 200)]


def test_concat_three_files(tmp_path):
    a = _write(tmp_path / "a.vcf", _PART_A)
    b = _write(tmp_path / "b.vcf", _PART_B)
    c = _write(tmp_path / "c.vcf", _HEADER + "chr2\t9000\t.\tT\tA\t60\tPASS\t.\tGT\t0/1\n")
    out = tmp_path / "out.vcf"
    n = bcftools.concat([str(a), str(b), str(c)], str(out))
    assert n == 4


def test_concat_one_input_raises(tmp_path):
    a = _write(tmp_path / "a.vcf", _PART_A)
    out = tmp_path / "out.vcf"
    with pytest.raises((IOError, OSError, ValueError)):
        bcftools.concat([str(a)], str(out))


def test_concat_incompatible_samples_raises(tmp_path):
    a = _write(tmp_path / "a.vcf", _PART_A)
    diff_samples = """##fileformat=VCFv4.3
##contig=<ID=chr1,length=10000>
##contig=<ID=chr2,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA99999
chr1\t900\t.\tA\tG\t30\tPASS\t.\tGT\t0/1
"""
    b = _write(tmp_path / "b.vcf", diff_samples)
    out = tmp_path / "out.vcf"
    with pytest.raises((IOError, OSError)):
        bcftools.concat([str(a), str(b)], str(out))


def test_concat_incompatible_contigs_raises(tmp_path):
    a = _write(tmp_path / "a.vcf", _PART_A)
    diff_contigs = """##fileformat=VCFv4.3
##contig=<ID=chrX,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878
chrX\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0/1
"""
    b = _write(tmp_path / "b.vcf", diff_contigs)
    out = tmp_path / "out.vcf"
    with pytest.raises((IOError, OSError)):
        bcftools.concat([str(a), str(b)], str(out))


def test_concat_output_bgzf(tmp_path):
    a = _write(tmp_path / "a.vcf", _PART_A)
    b = _write(tmp_path / "b.vcf", _PART_B)
    out = tmp_path / "out.vcf.gz"
    n = bcftools.concat([str(a), str(b)], str(out), output_type="z")
    assert n == 3
    assert out.read_bytes()[:2] == b"\x1f\x8b"


def test_concat_output_bcf(tmp_path):
    a = _write(tmp_path / "a.vcf", _PART_A)
    b = _write(tmp_path / "b.vcf", _PART_B)
    out = tmp_path / "out.bcf"
    n = bcftools.concat([str(a), str(b)], str(out), output_type="b")
    assert n == 3
    with rubam.VariantFile(str(out), "r") as f:
        assert sum(1 for _ in f) == 3
