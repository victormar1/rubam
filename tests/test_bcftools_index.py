"""Smoke + functional tests for rubam.tools.bcftools.index."""
import shutil
import subprocess
from pathlib import Path

import pytest

import rubam
from rubam.tools import bcftools
from _wsl_probe import wsl_usable


_TINY_VCF = """##fileformat=VCFv4.3
##contig=<ID=chr1,length=10000>
##contig=<ID=chr2,length=10000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30\tPASS\t.
chr1\t500\t.\tC\tT\t40\tPASS\t.
chr1\t1500\t.\tG\tA\t30\tPASS\t.
chr2\t200\t.\tA\tT\t30\tPASS\t.
"""


def _to_wsl(p: Path) -> str:
    s = str(p).replace("\\", "/")
    if len(s) > 1 and s[1] == ":":
        s = "/mnt/" + s[0].lower() + s[2:]
    return s


@pytest.fixture
def vcf_gz(tmp_path):
    """A bgzipped VCF without an index, ready for `bcftools index`."""
    plain = tmp_path / "tiny.vcf"
    plain.write_text(_TINY_VCF)
    gz = tmp_path / "tiny.vcf.gz"
    if not wsl_usable():
        pytest.skip("WSL Ubuntu (bgzip) not usable")
    subprocess.check_call(
        ["wsl", "-d", "Ubuntu", "bash", "-lc",
         f"bgzip -c {_to_wsl(plain)} > {_to_wsl(gz)}"]
    )
    return gz


def test_index_tbi_default(vcf_gz):
    out = bcftools.index(str(vcf_gz))
    assert out.endswith(".tbi")
    assert Path(out).exists()


def test_index_csi(vcf_gz):
    out = bcftools.index(str(vcf_gz), csi=True)
    assert out.endswith(".csi")
    assert Path(out).exists()


def test_index_then_fetch(vcf_gz):
    # Build TBI, then verify that VariantFile.fetch works.
    bcftools.index(str(vcf_gz))
    with rubam.VariantFile(str(vcf_gz), "r") as f:
        recs = list(f.fetch("chr1", 1, 1000))
    positions = sorted(r.position for r in recs)
    assert positions == [100, 500]


def test_index_existing_no_force_fails(vcf_gz):
    bcftools.index(str(vcf_gz))
    with pytest.raises((IOError, OSError, ValueError)):
        bcftools.index(str(vcf_gz))  # force=False -> error


def test_index_force_overwrites(vcf_gz):
    bcftools.index(str(vcf_gz))
    out = bcftools.index(str(vcf_gz), force=True)
    assert Path(out).exists()


def test_index_plain_vcf_rejected(tmp_path):
    plain = tmp_path / "tiny.vcf"
    plain.write_text(_TINY_VCF)
    with pytest.raises((IOError, OSError, ValueError)):
        bcftools.index(str(plain))


def test_index_bcf_tbi_rejected(tmp_path):
    # Use rubam to write a BCF, then attempt --tbi (should reject).
    src = tmp_path / "tiny.vcf"
    src.write_text(_TINY_VCF)
    bcf_path = tmp_path / "tiny.bcf"
    with rubam.VariantFile(str(src), "r") as r:
        with rubam.VariantFile(str(bcf_path), "wb", header=r.header) as w:
            for rec in r:
                w.write(rec)
    # csi=False (default TBI) on BCF must be rejected
    with pytest.raises((IOError, OSError, ValueError)):
        bcftools.index(str(bcf_path), csi=False)


def test_index_bcf_csi_works(tmp_path):
    src = tmp_path / "tiny.vcf"
    src.write_text(_TINY_VCF)
    bcf_path = tmp_path / "tiny.bcf"
    with rubam.VariantFile(str(src), "r") as r:
        with rubam.VariantFile(str(bcf_path), "wb", header=r.header) as w:
            for rec in r:
                w.write(rec)
    out = bcftools.index(str(bcf_path), csi=True)
    assert out.endswith(".csi")
    assert Path(out).exists()
