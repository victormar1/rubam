"""Windows-native stress tests for rubam.

Goal: prove rubam handles realistic Windows path/filesystem edge cases that a
clinical/hospital Windows user would actually hit. These tests are additive
(skipped on non-Windows) and exercise the *open path* via the real public API
(rubam.AlignmentFile / rubam.get_depths) — not just os.stat.

Each test copies the canonical fixture (tests/fixtures/smoke.bam + .bai) into a
short-lived TemporaryDirectory at a "weird" path, opens it, asserts something
non-trivial, then cleans up.

Source fixture:
    tests/fixtures/smoke.bam     — single-contig synthetic BAM (chr1, len=1000)
    tests/fixtures/smoke.bam.bai — index for the above

The get_depths API is 1-based inclusive (validated 2026-05-13): start>=1.
"""

from __future__ import annotations

import os
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

import rubam

WINDOWS_ONLY = pytest.mark.skipif(
    sys.platform != "win32",
    reason="Windows-native path stress test; irrelevant off Windows",
)

FIXTURE_BAM = Path(__file__).parent / "fixtures" / "smoke.bam"
FIXTURE_BAI = Path(__file__).parent / "fixtures" / "smoke.bam.bai"


def _copy_fixture_to(dest_dir: Path, filename: str = "smoke.bam") -> Path:
    """Copy smoke.bam + smoke.bam.bai into dest_dir under `filename` and return
    the absolute path to the BAM."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    bam_dst = dest_dir / filename
    bai_dst = dest_dir / (filename + ".bai")
    shutil.copy2(FIXTURE_BAM, bam_dst)
    shutil.copy2(FIXTURE_BAI, bai_dst)
    return bam_dst


def _assert_smoke_bam_works(bam_path: str | os.PathLike) -> None:
    """Open via AlignmentFile + run get_depths. Asserts the BAM is the smoke
    fixture (chr1, length 1000, non-zero coverage in [1,100])."""
    path_str = str(bam_path)
    with rubam.AlignmentFile(path_str, "rb") as bam:
        assert bam.is_open
        assert "chr1" in bam.references, f"references={bam.references}"
        assert bam.lengths[0] == 1000
    # Exercise the parallel-Rust read path too:
    positions, depths = rubam.get_depths(path_str, "chr1", 1, 100)
    assert len(positions) == 100
    assert len(depths) == 100
    assert max(depths) > 0, "expected non-zero coverage somewhere in [1,100]"


# ---------------------------------------------------------------------------
# a) Path containing spaces
# ---------------------------------------------------------------------------
@WINDOWS_ONLY
def test_path_with_spaces():
    """BAM at a path containing spaces — e.g. 'C:\\Users\\X\\Path With Space\\test.bam'."""
    with tempfile.TemporaryDirectory(prefix="rubam_stress_") as td:
        weird_dir = Path(td) / "Path With Space"
        bam = _copy_fixture_to(weird_dir, "test.bam")
        assert " " in str(bam)
        _assert_smoke_bam_works(bam)


# ---------------------------------------------------------------------------
# b) Path containing accents / non-ASCII characters
# ---------------------------------------------------------------------------
@WINDOWS_ONLY
def test_path_with_accents():
    """BAM at a path containing non-ASCII chars (é, ü, 中). NTFS is UTF-16 so
    this should work, but std::fs canonicalize-on-open behavior is worth a check."""
    with tempfile.TemporaryDirectory(prefix="rubam_stress_") as td:
        weird_dir = Path(td) / "é_ü_中"
        bam = _copy_fixture_to(weird_dir, "test.bam")
        assert any(ord(c) > 127 for c in str(bam))
        _assert_smoke_bam_works(bam)


# ---------------------------------------------------------------------------
# c) Forward-slash vs back-slash equivalence
# ---------------------------------------------------------------------------
@WINDOWS_ONLY
def test_forward_vs_backslash_equivalence():
    """A Windows-aware library must accept both 'D:/x/y.bam' and 'D:\\x\\y.bam'
    for the same underlying file."""
    with tempfile.TemporaryDirectory(prefix="rubam_stress_") as td:
        bam = _copy_fixture_to(Path(td), "test.bam")
        as_back = str(bam).replace("/", "\\")
        as_fwd = str(bam).replace("\\", "/")
        # Both forms should open identically.
        _assert_smoke_bam_works(as_back)
        _assert_smoke_bam_works(as_fwd)


# ---------------------------------------------------------------------------
# d) UNC-style / Win32 long-path prefix (\\?\)
# ---------------------------------------------------------------------------
@WINDOWS_ONLY
def test_unc_longpath_prefix():
    """The '\\\\?\\' prefix is the Win32 normalized-path syntax used to bypass
    MAX_PATH and reserved-name parsing. Hospital deployments routinely hit
    paths long enough to require it. If rubam can't open via this prefix we
    want to know."""
    with tempfile.TemporaryDirectory(prefix="rubam_stress_") as td:
        bam = _copy_fixture_to(Path(td), "smoke.bam")
        # Build the verbatim form: \\?\<absolute-drive-path-with-backslashes>
        abs_back = str(bam.resolve()).replace("/", "\\")
        if abs_back.startswith("\\\\?\\"):
            unc_path = abs_back
        else:
            unc_path = "\\\\?\\" + abs_back
        _assert_smoke_bam_works(unc_path)


# ---------------------------------------------------------------------------
# e) Relative path resolution (./smoke.bam)
# ---------------------------------------------------------------------------
@WINDOWS_ONLY
def test_relative_path_resolution():
    """rubam should resolve a relative path against the current working
    directory, like every other BAM toolkit.

    NOTE on Windows tempdir cleanup: if the cwd is *inside* the TemporaryDirectory
    when the context exits, rmtree fails with WinError 32 because the cwd itself
    is "in use". We therefore restore cwd manually BEFORE removing the dir,
    then delete it with shutil.rmtree(ignore_errors=True) — the cross-version
    equivalent of TemporaryDirectory's 3.10+ ignore_cleanup_errors=True, which
    is unavailable on the 3.8/3.9 interpreters we still support.
    """
    saved_cwd = os.getcwd()
    td = tempfile.mkdtemp(prefix="rubam_stress_")
    try:
        td_path = Path(td)
        _copy_fixture_to(td_path, "smoke.bam")
        os.chdir(td_path)
        _assert_smoke_bam_works("./smoke.bam")
        _assert_smoke_bam_works("smoke.bam")
    finally:
        os.chdir(saved_cwd)
        shutil.rmtree(td, ignore_errors=True)


# ---------------------------------------------------------------------------
# f) pathlib.Path object accepted directly
# ---------------------------------------------------------------------------
@WINDOWS_ONLY
def test_pathlib_path_object():
    """**v0.3.2** (Wave 11): both `AlignmentFile(path)` AND the depth APIs
    (`get_depths`, `get_depths_numpy`, `depth_chunks`) accept `pathlib.Path`
    natively. The PyO3 bindings take `std::path::PathBuf`, which derives
    `FromPyObject` over `os.PathLike`."""
    with tempfile.TemporaryDirectory(prefix="rubam_stress_") as td:
        bam = _copy_fixture_to(Path(td), "smoke.bam")
        assert isinstance(bam, Path)
        # AlignmentFile with a real Path object:
        with rubam.AlignmentFile(bam, "rb") as f:
            assert f.is_open
            assert "chr1" in f.references
        # get_depths with pathlib.Path:
        positions, depths = rubam.get_depths(bam, "chr1", 1, 50)
        assert len(positions) == 50 and max(depths) >= 0
        # get_depths_numpy with pathlib.Path:
        positions_np, depths_np = rubam.get_depths_numpy(bam, "chr1", 1, 50)
        assert positions_np.shape == (50,) and depths_np.shape == (50,)
        # depth_chunks with pathlib.Path:
        chunks = list(rubam.depth_chunks(bam, "chr1", 1, 100, chunk_size=30))
        assert len(chunks) >= 1
        for pos_c, dep_c in chunks:
            assert pos_c.shape == dep_c.shape


# ---------------------------------------------------------------------------
# g) Open + close + reopen same file in the same process
# ---------------------------------------------------------------------------
@WINDOWS_ONLY
def test_open_close_reopen_same_file():
    """Windows is famously strict about file-handle lifecycles. Reopening the
    same BAM in the same process after a clean close must work — clinical
    pipelines do this routinely (per-sample loops)."""
    with tempfile.TemporaryDirectory(prefix="rubam_stress_") as td:
        bam = _copy_fixture_to(Path(td), "smoke.bam")
        path_str = str(bam)

        # First open
        f1 = rubam.AlignmentFile(path_str, "rb")
        assert f1.is_open
        refs1 = f1.references
        f1.close()
        assert not f1.is_open

        # Reopen in the same process
        f2 = rubam.AlignmentFile(path_str, "rb")
        assert f2.is_open
        refs2 = f2.references
        f2.close()
        assert not f2.is_open

        assert refs1 == refs2

        # Bonus: third reopen through the context manager API.
        with rubam.AlignmentFile(path_str, "rb") as f3:
            assert f3.is_open
            assert f3.references == refs1
