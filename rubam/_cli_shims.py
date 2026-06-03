"""pysam-compatible subprocess shims for `bcftools` and `samtools`.

`pysam.bcftools` and `pysam.samtools` are themselves subprocess wrappers
around the system `bcftools` / `samtools` binaries — the variant-caller
(`bcftools mpileup` / `bcftools call`), region-extraction (`samtools view`
with complex flag combinations) and many other utilities are NOT
re-implemented in pure Python by pysam; pysam shells out and returns the
stdout. We mirror that contract here.

Resolution order for each call:

1. **System binary** (`shutil.which("bcftools")` or `which("samtools")`)
   — first preference when available, matches pysam byte-for-byte on
   stdout/stderr.
2. **Bundled rubam companion binary** (`rubam-bcftools` or
   `rubam-samtools`) — fallback on Windows-only hosts that do not ship a
   system `bcftools` / `samtools`. Covers the documented subset of
   subcommands; raises a clear error on the rest.
3. **`NotImplementedError`** with a pointer to install instructions if
   neither is found.

Usage matches `pysam.bcftools.<subcmd>(...)`:

    >>> rubam.bcftools("mpileup", "-f", ref, bam)            # noqa: F821
    b'...'  # raw stdout bytes
    >>> rubam.bcftools("call", "-c", "-v", vcf_in)            # noqa: F821
    b'...'
    >>> rubam.samtools("view", "-c", "-F", "256", bam)        # noqa: F821
    b'347291\\n'

The function returns the subprocess' raw stdout bytes (not text), to
match pysam's `pysam.bcftools(*argv)` convention. Set the keyword
`catch_stdout=False` to stream directly to the parent's stdout (useful
for very large outputs).
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional


# Subcommands the rubam-bcftools companion can handle.
_RUBAM_BCFTOOLS_OK = {"view", "query", "sort", "index", "head", "concat", "stats", "norm"}
# Subcommands the rubam-samtools companion can handle.
_RUBAM_SAMTOOLS_OK = {"view", "index", "flagstat", "coverage", "idxstats", "head", "sort", "merge", "faidx", "calmd"}


def _bundled_binary(name: str) -> Optional[Path]:
    """Locate `rubam-bcftools` / `rubam-samtools` shipped with the wheel,
    or in the dev source tree's target/release/ directory."""
    # 1. Same directory as the rubam Python package (wheel-bundled).
    pkg_dir = Path(__file__).resolve().parent
    candidate = pkg_dir / (f"{name}.exe" if sys.platform == "win32" else name)
    if candidate.exists():
        return candidate
    # 2. Source tree `target/release/<name>`.
    src_target = pkg_dir.parent / "target" / "release" / (f"{name}.exe" if sys.platform == "win32" else name)
    if src_target.exists():
        return src_target
    # 3. PATH.
    p = shutil.which(name)
    return Path(p) if p else None


def _run(argv: list[str], *, catch_stdout: bool = True, check: bool = True) -> bytes:
    """Run a subprocess; return its raw stdout (or empty bytes if
    catch_stdout=False). Raises subprocess.CalledProcessError on
    non-zero exit when check=True."""
    if catch_stdout:
        proc = subprocess.run(argv, capture_output=True, check=False)
        if check and proc.returncode != 0:
            err = proc.stderr.decode("utf-8", errors="replace")
            raise subprocess.CalledProcessError(
                proc.returncode, argv, output=proc.stdout, stderr=proc.stderr,
            )
        return proc.stdout
    else:
        proc = subprocess.run(argv, check=False)
        if check and proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, argv)
        return b""


def _resolve_backend(system_name: str, companion_name: str,
                     companion_supports: set[str], subcmd: str) -> tuple[str, str]:
    """Return (binary_path, kind) where kind is 'system' or 'rubam'."""
    sys_bin = shutil.which(system_name)
    if sys_bin:
        return sys_bin, "system"
    bundled = _bundled_binary(companion_name)
    if bundled is not None and subcmd in companion_supports:
        return str(bundled), "rubam"
    raise NotImplementedError(
        f"{system_name}: neither the system `{system_name}` binary is on PATH "
        f"nor does the bundled `{companion_name}` companion implement `{subcmd}` "
        f"(supported by rubam: {sorted(companion_supports)}). "
        f"Install `{system_name}` from your package manager (Linux: "
        f"`apt install {system_name}` / `brew install {system_name}` / "
        f"`conda install -c bioconda {system_name}`; on Windows use WSL or "
        f"a static {system_name} build) or open a rubam issue if the subcommand "
        f"belongs in the companion's scope."
    )


def bcftools(*argv: str, catch_stdout: bool = True, check: bool = True) -> bytes:
    """`pysam.bcftools(*argv)` drop-in.

    First argument is the subcommand (`mpileup`, `call`, `view`, `query`,
    `sort`, `index`, `head`, `norm`, `stats`, ...). Remaining arguments
    are passed through as-is. Returns the subprocess' raw stdout bytes
    by default.

    For subcommands the rubam companion does not implement (notably
    `mpileup` and `call` — variant calling is out of the rubam in-scope),
    the system `bcftools` binary is used. If neither is present, a
    `NotImplementedError` is raised with install instructions.
    """
    if not argv:
        raise ValueError("rubam.bcftools requires a subcommand as the first argument")
    subcmd = argv[0]
    binary, _kind = _resolve_backend("bcftools", "rubam-bcftools", _RUBAM_BCFTOOLS_OK, subcmd)
    cmd = [binary, *argv]
    return _run(cmd, catch_stdout=catch_stdout, check=check)


def samtools(*argv: str, catch_stdout: bool = True, check: bool = True) -> bytes:
    """`pysam.samtools(*argv)` drop-in.

    First argument is the subcommand (`view`, `index`, `flagstat`,
    `idxstats`, `sort`, `merge`, `faidx`, `calmd`, `depth`, ...).
    Returns the subprocess' raw stdout bytes by default.

    For subcommands the rubam companion does not implement, the system
    `samtools` binary is used. The bundled `rubam-samtools` companion
    handles most of the inspection surface (see
    `docs/samtools_compatibility.md`).
    """
    if not argv:
        raise ValueError("rubam.samtools requires a subcommand as the first argument")
    subcmd = argv[0]
    # samtools depth is wired as a standalone binary `rubam-depth`, not as a
    # subcommand of rubam-samtools — bridge that here.
    if subcmd == "depth":
        sys_bin = shutil.which("samtools")
        if sys_bin:
            return _run([sys_bin, *argv], catch_stdout=catch_stdout, check=check)
        depth_bin = _bundled_binary("rubam-depth")
        if depth_bin is not None and len(argv) >= 4:
            # rubam-depth expects: <bam> <chrom> <start> <end> [flags...]
            # samtools depth expects: [flags...] [bam ...]
            # Defer to system samtools' calling convention if possible;
            # otherwise raise a clear error and point at rubam.depth().
            raise NotImplementedError(
                "rubam.samtools('depth', ...): the bundled rubam-depth binary has a "
                "different calling convention than `samtools depth`. Use "
                "`rubam.depth(bam, region=...)` for the pysam-compatible Python "
                "API, or install system `samtools` for byte-for-byte CLI parity."
            )
    binary, _kind = _resolve_backend("samtools", "rubam-samtools", _RUBAM_SAMTOOLS_OK, subcmd)
    cmd = [binary, *argv]
    return _run(cmd, catch_stdout=catch_stdout, check=check)


# Backwards-compat shorthand: pysam exposes `pysam.bcftools.mpileup(...)`
# via attribute-style access. Mirror that.
class _SubcommandDispatcher:
    """`rubam.bcftools.mpileup(*args)` shortcut — equivalent to
    `rubam.bcftools("mpileup", *args)`."""
    def __init__(self, fn):
        self._fn = fn

    def __call__(self, *argv, **kwargs):
        return self._fn(*argv, **kwargs)

    def __getattr__(self, subcmd: str):
        def _bound(*args, **kwargs):
            return self._fn(subcmd, *args, **kwargs)
        _bound.__name__ = f"{self._fn.__name__}.{subcmd}"
        _bound.__doc__ = f"Shortcut for `{self._fn.__name__}({subcmd!r}, *args)`."
        return _bound


bcftools = _SubcommandDispatcher(bcftools)
samtools = _SubcommandDispatcher(samtools)
