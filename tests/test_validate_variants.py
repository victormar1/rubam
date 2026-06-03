"""Smoke test for the cross-tool VCF validator.

Runs bench.validate_variants inside WSL against the synthetic 3-sample
100-record fixture. Skipped if:
  - WSL is not available on this machine, or
  - tests/data/validation_3sample_100rec.vcf.gz does not exist
    (run scripts/build_validation_vcf.sh first).
"""
import shutil
import subprocess
from pathlib import Path

import pytest

from _wsl_probe import wsl_usable


VCF = Path("tests/data/validation_3sample_100rec.vcf.gz")


def _to_wsl(p: Path) -> str:
    """Convert a Windows path to its WSL /mnt/... equivalent."""
    s = str(p).replace("\\", "/")
    if len(s) > 1 and s[1] == ":":
        s = "/mnt/" + s[0].lower() + s[2:]
    return s


@pytest.mark.skipif(not VCF.exists(), reason="run scripts/build_validation_vcf.sh first")
@pytest.mark.skipif(not wsl_usable(), reason="WSL Ubuntu (system pysam) not usable")
def test_rubam_vs_pysam_synth_3sample_100rec():
    rubam_dir = _to_wsl(Path.cwd().resolve())
    vcf_wsl = _to_wsl(VCF)
    cmd = (
        f"cd {rubam_dir} && "
        f"source ~/.rubam-venv/bin/activate && "
        f"python -m bench.validate_variants {vcf_wsl}"
    )
    result = subprocess.run(
        ["wsl", "-d", "Ubuntu", "bash", "-lc", cmd],
        capture_output=True,
        text=True,
    )
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)
    assert result.returncode == 0, f"validation failed:\n{result.stderr}"
    assert "100/100 records match" in result.stdout or "100/100" in result.stdout
