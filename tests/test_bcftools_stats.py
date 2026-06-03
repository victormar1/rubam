"""Smoke + functional tests for rubam.tools.bcftools.stats."""
import pytest
from rubam.tools import bcftools


_VCF_FOR_STATS = """##fileformat=VCFv4.3
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\tNA12891
chr1\t100\t.\tA\tG\t30\t.\t.\tGT\t0/0\t0/1
chr1\t200\t.\tC\tT\t30\t.\t.\tGT\t0/1\t1/1
chr1\t300\t.\tG\tA\t30\t.\t.\tGT\t1/1\t0/0
chr1\t400\t.\tA\tT\t30\t.\t.\tGT\t0/1\t0/0
chr1\t500\t.\tACGT\tA\t30\t.\t.\tGT\t0/1\t0/0
chr1\t600\t.\tAA\tCC\t30\t.\t.\tGT\t0/0\t0/1
"""


@pytest.fixture
def src_vcf(tmp_path):
    p = tmp_path / "src.vcf"
    p.write_text(_VCF_FOR_STATS)
    return p


def test_stats_total_records(src_vcf):
    s = bcftools.stats(str(src_vcf))
    assert s["total_records"] == 6


def test_stats_classification(src_vcf):
    s = bcftools.stats(str(src_vcf))
    # SNPs: positions 100 (A>G), 200 (C>T), 300 (G>A), 400 (A>T) = 4
    # Indel: position 500 (ACGT>A) = 1
    # MNP: position 600 (AA>CC) = 1
    assert s["snps"] == 4
    assert s["indels"] == 1
    assert s["mnps"] == 1
    assert s["complex"] == 0


def test_stats_ts_tv(src_vcf):
    s = bcftools.stats(str(src_vcf))
    # Transitions: A>G, C>T, G>A = 3
    # Transversions: A>T = 1
    assert s["transitions"] == 3
    assert s["transversions"] == 1
    assert abs(s["ts_tv_ratio"] - 3.0) < 1e-9


def test_stats_sample_counts(src_vcf):
    s = bcftools.stats(str(src_vcf))
    samples = s["samples"]
    assert set(samples.keys()) == {"NA12878", "NA12891"}
    # NA12878 GTs: 0/0 0/1 1/1 0/1 0/1 0/0
    #    = 2 hom_ref, 3 het, 1 hom_alt
    assert samples["NA12878"]["hom_ref"] == 2
    assert samples["NA12878"]["het"] == 3
    assert samples["NA12878"]["hom_alt"] == 1
    assert samples["NA12878"]["missing"] == 0
    # NA12891 GTs: 0/1 1/1 0/0 0/0 0/0 0/1
    #    = 3 hom_ref, 2 het, 1 hom_alt
    assert samples["NA12891"]["hom_ref"] == 3
    assert samples["NA12891"]["het"] == 2
    assert samples["NA12891"]["hom_alt"] == 1
    assert samples["NA12891"]["missing"] == 0


def test_stats_missing_genotype(tmp_path):
    p = tmp_path / "missing.vcf"
    p.write_text("""##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
chr1\t100\t.\tA\tG\t30\t.\t.\tGT\t./.
""")
    s = bcftools.stats(str(p))
    assert s["samples"]["S"]["missing"] == 1
    assert s["samples"]["S"]["hom_ref"] == 0


def test_stats_empty_vcf(tmp_path):
    p = tmp_path / "empty.vcf"
    p.write_text("""##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT
""")
    s = bcftools.stats(str(p))
    assert s["total_records"] == 0
    assert s["samples"] == {}
    assert s["ts_tv_ratio"] == 0.0
