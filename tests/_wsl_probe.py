"""Shared probe for tests that cross-validate against system bioinformatics
tools (samtools / bcftools / bgzip) run through WSL on Windows.

A bare `shutil.which("wsl")` check is NOT enough: GitHub's `windows-latest`
runner ships `wsl.exe` but has no Linux distro installed, so `wsl -d Ubuntu …`
fails at runtime. These cross-validation tests must SKIP (not error/fail) when
a usable WSL Ubuntu environment is absent, so they only run on a developer
machine that actually has the tools.
"""
from __future__ import annotations

import functools
import shutil
import subprocess


@functools.lru_cache(maxsize=1)
def wsl_usable() -> bool:
    """True only if WSL exists AND an `Ubuntu` distro can run a shell."""
    if shutil.which("wsl") is None:
        return False
    try:
        proc = subprocess.run(
            ["wsl", "-d", "Ubuntu", "bash", "-lc", "true"],
            capture_output=True,
            timeout=30,
        )
        return proc.returncode == 0
    except Exception:
        return False
