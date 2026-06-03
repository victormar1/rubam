"""Smoke + functional tests for rubam.tools.bcftools.query."""
import pytest
from rubam.tools import bcftools

_VCF = """##fileformat=VCFv4.3
##contig=<ID=chr1,length=10000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\tNA12891
chr1\t100\trs1\tA\tG\t30\tPASS\tDP=20;AF=0.5\tGT:DP:AD\t0/1:30:15,15\t1/1:25:0,25
chr1\t500\t.\tC\tT\t40\tPASS\tDP=15\tGT:DP:AD\t0/0:15:15,0\t0/1:20:10,10
"""

@pytest.fixture
def src_vcf(tmp_path):
    p = tmp_path / "src.vcf"
    p.write_text(_VCF)
    return p


def test_query_chrom_pos(src_vcf, tmp_path):
    out = tmp_path / "out.txt"
    n = bcftools.query(str(src_vcf), r"%CHROM\t%POS\n", output=str(out))
    assert n == 2
    assert out.read_text() == "chr1\t100\nchr1\t500\n"


def test_query_ref_alt(src_vcf, tmp_path):
    out = tmp_path / "out.txt"
    n = bcftools.query(str(src_vcf), r"%CHROM\t%POS\t%REF\t%ALT\n", output=str(out))
    assert n == 2
    lines = out.read_text().strip().split("\n")
    assert lines[0] == "chr1\t100\tA\tG"
    assert lines[1] == "chr1\t500\tC\tT"


def test_query_id_missing(src_vcf, tmp_path):
    out = tmp_path / "out.txt"
    n = bcftools.query(str(src_vcf), r"%ID\n", output=str(out))
    assert n == 2
    lines = out.read_text().strip().split("\n")
    assert lines[0] == "rs1"
    assert lines[1] == "."  # missing ID rendered as .


def test_query_info_field(src_vcf, tmp_path):
    out = tmp_path / "out.txt"
    n = bcftools.query(str(src_vcf), r"%CHROM\t%POS\t%INFO/DP\n", output=str(out))
    assert n == 2
    lines = out.read_text().strip().split("\n")
    assert lines[0].endswith("\t20")
    assert lines[1].endswith("\t15")


def test_query_qual_and_filter(src_vcf, tmp_path):
    out = tmp_path / "out.txt"
    n = bcftools.query(str(src_vcf), r"%QUAL\t%FILTER\n", output=str(out))
    assert n == 2
    text = out.read_text()
    assert "30" in text and "PASS" in text


def test_query_per_sample_gt(src_vcf, tmp_path):
    out = tmp_path / "out.txt"
    n = bcftools.query(
        str(src_vcf),
        r"%CHROM\t%POS[\t%SAMPLE=%GT]\n",
        output=str(out),
    )
    assert n == 2
    line0 = out.read_text().split("\n")[0]
    # Order of samples comes from the header
    assert "NA12878=0/1" in line0
    assert "NA12891=1/1" in line0


def test_query_per_sample_dp(src_vcf, tmp_path):
    out = tmp_path / "out.txt"
    n = bcftools.query(
        str(src_vcf),
        r"[%SAMPLE %DP\n]",
        output=str(out),
    )
    assert n == 2
    text = out.read_text()
    # Two records × 2 samples = 4 lines
    lines = text.strip().split("\n")
    assert len(lines) == 4
    assert "NA12878 30" in lines or "NA12878 30" in text
    assert "NA12891 25" in lines or "NA12891 25" in text


def test_query_array_field_ad(src_vcf, tmp_path):
    out = tmp_path / "out.txt"
    n = bcftools.query(
        str(src_vcf),
        r"[%AD\n]",
        output=str(out),
    )
    assert n == 2
    text = out.read_text()
    assert "15,15" in text  # array joined with ,
    assert "0,25" in text


def test_query_unknown_info_key_renders_dot(src_vcf, tmp_path):
    out = tmp_path / "out.txt"
    n = bcftools.query(str(src_vcf), r"%INFO/MISSING\n", output=str(out))
    assert n == 2
    lines = out.read_text().strip().split("\n")
    assert lines == [".", "."]


def test_query_empty_format_raises(src_vcf, tmp_path):
    with pytest.raises(ValueError):
        bcftools.query(str(src_vcf), "", output=str(tmp_path / "x.txt"))
