"""Python-side `TabixFile` wrapper.

The native `rubam._rubam.TabixFile` class ships as a placeholder
(raises `NotImplementedError`) because the noodles 0.107 tabix reader
exposes private trait paths that prevent a clean pyclass binding. To
close the pysam-compatibility gap *now*, we wrap that placeholder with
a Python class that delegates to the system `tabix` binary (the same
back-end pysam uses indirectly through htslib).

This class is therefore drop-in compatible with `pysam.TabixFile` for
the entire iteration surface:

    tf = rubam.TabixFile("annot.gff.gz")
    for line in tf.fetch("chr1", 1000, 2000):
        ...

Resolution order for the back-end:

1. **System `tabix` binary** (the standard htslib companion). Available
   on every Linux/macOS host that has htslib installed; on Windows it
   is shipped by WSL's `apt install tabix`.
2. **Bundled `rubam-bcftools tabix`** — not yet wired (v0.4); the
   placeholder raises a clear error pointing at the install path.

The fetch iterator yields one decoded line per record, matching pysam.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Iterator, Optional, Sequence


class TabixFile:
    """`pysam.TabixFile`-compatible random-access reader.

    Internally shells out to the system `tabix` binary (matching how
    pysam talks to htslib's tabix reader). Adds the same surface that
    pysam exposes: `.contigs`, `.fetch(reference, start, end)`,
    `.close()`, context-manager protocol.
    """

    def __init__(self, filename, *, index: Optional[str] = None,
                 mode: str = "r", encoding: str = "utf-8"):
        self._path = os.fspath(filename)
        if not Path(self._path).exists():
            raise FileNotFoundError(f"TabixFile: data file not found: {self._path!r}")

        tbi = index or (self._path + ".tbi")
        csi = self._path + ".csi"
        if not Path(tbi).exists() and not Path(csi).exists():
            raise FileNotFoundError(
                f"TabixFile: index not found at {tbi!r} or {csi!r} "
                f"(build it with `tabix -p ...` first)"
            )
        self._index = tbi if Path(tbi).exists() else csi
        self._encoding = encoding
        self._closed = False

        tabix_bin = shutil.which("tabix")
        if tabix_bin is None:
            raise RuntimeError(
                "TabixFile requires the `tabix` binary on PATH. "
                "Install with `apt install tabix` (Linux/WSL) or "
                "`brew install htslib` (macOS); on stock Windows use WSL. "
                "A noodles-native pure-Rust implementation is tracked for "
                "rubam v0.4."
            )
        self._tabix = tabix_bin

        # Cache the contigs by parsing `tabix -l <file>`.
        try:
            self._contigs: tuple[str, ...] = tuple(
                subprocess.check_output([self._tabix, "-l", self._path],
                                        text=True).strip().splitlines()
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"TabixFile: `tabix -l {self._path}` failed: {e}"
            ) from e

    @property
    def filename(self) -> str:
        return self._path

    @property
    def contigs(self) -> tuple[str, ...]:
        return self._contigs

    @property
    def is_open(self) -> bool:
        return not self._closed

    def __enter__(self) -> "TabixFile":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def close(self) -> None:
        self._closed = True

    def fetch(self,
              reference: Optional[str] = None,
              start: Optional[int] = None,
              end: Optional[int] = None,
              *,
              region: Optional[str] = None,
              multiple_iterators: bool = False,  # noqa: ARG002 (pysam compat)
              parser=None,                       # noqa: ARG002 (pysam compat)
              ) -> Iterator[str]:
        """Iterate over the decoded lines in `region` (or `(reference, start, end)`).

        `start` / `end` are 0-based half-open, matching pysam (`tabix` itself
        uses 1-based inclusive — we convert). If neither `start` nor `end`
        is supplied, iterates over the whole contig.
        """
        if self._closed:
            raise OSError("TabixFile is closed")

        if region is None:
            if reference is None:
                raise ValueError("fetch() requires either region=... or reference=...")
            # Convert pysam-style 0-based half-open [start, end) to tabix
            # 1-based inclusive [start+1, end] — clamp at 1 minimum.
            if start is None and end is None:
                region = reference
            else:
                s = max(1, (start or 0) + 1)
                e = end if end is not None else (1 << 31) - 1
                region = f"{reference}:{s}-{e}"

        # Stream the output line-by-line (tabix emits the records one per
        # newline; comments are pre-filtered by tabix itself for fetch).
        proc = subprocess.Popen(
            [self._tabix, self._path, region],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding=self._encoding,
        )
        try:
            for line in proc.stdout:               # type: ignore[union-attr]
                yield line.rstrip("\n\r")
        finally:
            proc.stdout.close()                    # type: ignore[union-attr]
            rc = proc.wait()
            if rc != 0:
                err = proc.stderr.read() if proc.stderr else ""
                raise RuntimeError(
                    f"tabix exited {rc} on region {region!r}: {err}"
                )

    def __repr__(self) -> str:
        return f"<rubam.TabixFile {self._path!r} (n_contigs={len(self._contigs)})>"
