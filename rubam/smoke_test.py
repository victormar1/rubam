"""Minimal smoke test for an installed rubam wheel.

Run as: ``python -m rubam.smoke_test``

The point is to validate that the *wheel* itself is functional in a clean
environment, i.e. without the rubam source tree on ``sys.path``. It must
therefore rely exclusively on what the wheel ships:

  - the native extension (``rubam._rubam``)
  - the package data file ``tests/fixtures/smoke.bam`` (+ ``.bai``)
    accessed via :mod:`importlib.resources`.

Exit code is 0 on success, non-zero on any failure (assertion or
exception propagates).
"""

from __future__ import annotations

import sys
from importlib import resources
from pathlib import Path

import rubam


def _resolve_smoke_bam() -> Path:
    """Locate the bundled ``smoke.bam`` shipped inside the wheel.

    The fixture is included via the ``[tool.maturin]`` ``include`` directive
    (see ``pyproject.toml``). When the wheel is installed the file lives
    next to the ``rubam`` package, under ``tests/fixtures/smoke.bam``.

    We try ``importlib.resources`` first (the canonical, wheel-friendly
    path), and fall back to a filesystem walk anchored on the package
    install location so that the smoke test stays useful even if the
    packaging mechanics change.
    """

    # Preferred: importlib.resources on the rubam package itself.
    pkg_root = resources.files("rubam")
    candidate = pkg_root.joinpath("..", "tests", "fixtures", "smoke.bam")
    try:
        as_path = Path(str(candidate)).resolve()
        if as_path.exists():
            return as_path
    except (FileNotFoundError, NotImplementedError):
        pass

    # Fallback: walk a few common roots adjacent to the package.
    pkg_dir = Path(str(pkg_root)).resolve()
    for base in (pkg_dir.parent, pkg_dir.parent.parent):
        guess = base / "tests" / "fixtures" / "smoke.bam"
        if guess.exists():
            return guess

    raise FileNotFoundError(
        "smoke.bam fixture not found next to the installed rubam package. "
        "Check the [tool.maturin] include directive in pyproject.toml."
    )


def main() -> int:
    print(f"rubam version: {rubam.__version__}")

    bam_path = _resolve_smoke_bam()
    print(f"smoke bam: {bam_path}")

    # Synth BAM uses chr1 over [1, 1000]. Query a slice that's guaranteed
    # to have coverage given the synth params (5x over 1000 bp, 100 bp reads).
    positions, depths = rubam.get_depths(str(bam_path), "chr1", 100, 500)

    assert len(positions) == len(depths), "positions/depths length mismatch"
    assert len(positions) > 0, "expected at least one sampled position"
    assert any(d > 0 for d in depths), (
        "expected at least one position with non-zero depth; "
        f"got depths={depths[:20]}..."
    )

    print(f"sampled {len(positions)} positions, max depth = {max(depths)}")
    print("smoke test OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
