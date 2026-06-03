"""Acceptance tests for the v0.3.2 BAM write path (P0-1).

Closes the pysam-write gap that blocked downstream ports. Tests four patterns:

  1. ``AlignmentFile(path, "wb", template=other)`` — copies the source
     file's header into the output.
  2. ``AlignmentFile(path, "wb", header=hdr)`` — explicit header kwarg.
  3. Write-after-filter (the common pattern: read indexed, conditionally
     emit; no record modification).
  4. ``write()`` returns 0 on a read-mode file → must raise ``ValueError``.

Plus a roundtrip: write 50 records, index the output, reopen, count.

The output BAMs are also cross-validated with system ``samtools view -c``
under WSL when available. Skipped wholesale on non-Windows.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile

import pytest

import rubam


SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "smoke.bam")


def _has_smoke_bam() -> bool:
    return os.path.exists(SRC)


@pytest.fixture
def tmp_bam(tmp_path):
    return str(tmp_path / "out.bam")


@pytest.mark.skipif(not _has_smoke_bam(), reason="smoke.bam fixture missing")
def test_write_with_template(tmp_bam):
    bam_in = rubam.AlignmentFile(SRC, "rb")
    bam_out = rubam.AlignmentFile(tmp_bam, "wb", template=bam_in)
    n = 0
    for r in iter(bam_in):
        bam_out.write(r)
        n += 1
    bam_out.close()
    bam_in.close()
    assert n == 50, f"expected 50 records, wrote {n}"
    assert os.path.exists(tmp_bam) and os.path.getsize(tmp_bam) > 0


@pytest.mark.skipif(not _has_smoke_bam(), reason="smoke.bam fixture missing")
def test_write_with_explicit_header(tmp_bam):
    bam_in = rubam.AlignmentFile(SRC, "rb")
    bam_out = rubam.AlignmentFile(tmp_bam, "wb", header=bam_in.header)
    n = 0
    for r in iter(bam_in):
        bam_out.write(r)
        n += 1
    bam_out.close()
    bam_in.close()
    assert n == 50


@pytest.mark.skipif(not _has_smoke_bam(), reason="smoke.bam fixture missing")
def test_write_filter_pass_through(tmp_bam):
    """Filter-shape pattern — read indexed, conditionally pass-through."""
    bam_in = rubam.AlignmentFile(SRC, "rb")
    bam_out = rubam.AlignmentFile(tmp_bam, "wb", template=bam_in)
    kept = 0
    for r in iter(bam_in):
        if r.mapping_quality >= 60:
            bam_out.write(r)
            kept += 1
    bam_out.close()
    bam_in.close()
    # smoke.bam has all MAPQ=60, so all 50 pass through.
    assert kept == 50


@pytest.mark.skipif(not _has_smoke_bam(), reason="smoke.bam fixture missing")
def test_write_roundtrip(tmp_bam):
    """Write 50 records, index, reopen, count — must round-trip."""
    bam_in = rubam.AlignmentFile(SRC, "rb")
    bam_out = rubam.AlignmentFile(tmp_bam, "wb", template=bam_in)
    for r in iter(bam_in):
        bam_out.write(r)
    bam_out.close()
    bam_in.close()

    rubam.index(tmp_bam)
    assert os.path.exists(tmp_bam + ".bai")

    bam_re = rubam.AlignmentFile(tmp_bam, "rb")
    n = sum(1 for _ in iter(bam_re))
    bam_re.close()
    assert n == 50


def test_write_missing_template_and_header_errors(tmp_bam):
    with pytest.raises(ValueError, match=r"template|header"):
        rubam.AlignmentFile(tmp_bam, "wb")


@pytest.mark.skipif(not _has_smoke_bam(), reason="smoke.bam fixture missing")
def test_write_on_read_mode_errors():
    bam_in = rubam.AlignmentFile(SRC, "rb")
    first = next(iter(bam_in))
    with pytest.raises(ValueError, match=r"not opened for writing"):
        bam_in.write(first)
    bam_in.close()


@pytest.mark.skipif(
    sys.platform != "win32" or shutil.which("wsl") is None,
    reason="WSL not available — cross-validation against system samtools requires it",
)
def test_write_cross_validation_with_system_samtools(tmp_bam):
    """The gold-standard test: write a BAM with rubam, then count its records
    with system samtools (in WSL) and compare. Catches BGZF-EOF / header /
    record-encoding bugs that pure-Rust roundtrip might miss."""
    if not _has_smoke_bam():
        pytest.skip("smoke.bam missing")

    bam_in = rubam.AlignmentFile(SRC, "rb")
    bam_out = rubam.AlignmentFile(tmp_bam, "wb", template=bam_in)
    for r in iter(bam_in):
        bam_out.write(r)
    bam_out.close()
    bam_in.close()

    drive = tmp_bam[0].lower()
    wsl_path = "/mnt/" + drive + tmp_bam[2:].replace("\\", "/")
    proc = subprocess.run(
        ["wsl", "-d", "Ubuntu", "samtools", "view", "-c", wsl_path],
        capture_output=True, text=True, timeout=30,
    )
    assert proc.returncode == 0, f"samtools rejected the BAM: {proc.stderr}"
    assert proc.stdout.strip() == "50", (
        f"samtools view -c returned {proc.stdout!r}, expected '50'"
    )
