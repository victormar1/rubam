# Windows-native stress tests

**Last updated: 2026-05-14 (against rubam v0.3.2 / Python 3.11.10).**

Lightweight pytest suite at `tests/test_windows_native_stress.py` that
validates `rubam`'s behavior under the path/filesystem conditions a real
clinical/hospital Windows user encounters — not just the "pip install
succeeds" baseline.

Every case is auto-skipped on non-Windows, and every weird path is created
inside a `tempfile.TemporaryDirectory` so nothing is left on disk after the run.

## What is tested

All cases exercise the **real public API** — `rubam.AlignmentFile(...)`
*and* `rubam.get_depths(...)` — so the file is genuinely opened and
parsed, not just `os.stat`-ed.

| # | Case                                  | What it proves                                                                                                  |
|---|---------------------------------------|-----------------------------------------------------------------------------------------------------------------|
| a | Path containing **spaces**            | `C:\Users\X\Path With Space\test.bam` opens and yields the expected references / depths.                       |
| b | Path containing **accents / non-ASCII** | `<tmp>\é_ü_中\test.bam` opens (NTFS is UTF-16 — but Rust's `std::fs` canonicalization is worth a check).        |
| c | **Forward-slash vs back-slash** equivalence | The same file opens via `D:/...` and `D:\\...`.                                                                |
| d | **Win32 `\\?\` long-path prefix**     | Verbatim Win32 path form (used to bypass `MAX_PATH` in deep clinical pipelines) opens correctly.                |
| e | **Relative path** resolution          | `./smoke.bam` and `smoke.bam` resolve against the current working directory.                                    |
| f | **`pathlib.Path`** object accepted directly | API ergonomics parity with pysam (which accepts `Path` everywhere).                                            |
| g | **Open / close / reopen** in same process | Per-sample loops can reopen the same BAM after a clean `close()` — no Windows handle-lifecycle quirks.         |

## Source fixture

- `tests/fixtures/smoke.bam`     — single-contig synthetic BAM (`chr1`, length 1000)
- `tests/fixtures/smoke.bam.bai` — companion index

Copied into each test's `TemporaryDirectory` under the weird path; cleaned up
on test exit.

## How to run

```powershell
# from the repo root, inside your activated venv
python -m pytest tests/test_windows_native_stress.py -v
```

## Results — rubam v0.3.2 (Windows 11, Python 3.11.10)

```
7 passed
```

| # | Case                                   | Outcome  |
|---|----------------------------------------|----------|
| a | Path with spaces                       | PASS     |
| b | Path with accents / non-ASCII          | PASS     |
| c | Forward-slash vs back-slash            | PASS     |
| d | Win32 `\\?\` long-path prefix          | PASS     |
| e | Relative path resolution               | PASS     |
| f | `pathlib.Path` accepted directly       | PASS     |
| g | Open / close / reopen same file        | PASS     |

All seven cases — including the gnarly ones (accents, `\\?\` prefix,
forward/backslash mixing, native `pathlib.Path` acceptance) — pass
cleanly on the current build. PyO3 signatures now use `PathBuf`,
which derives `FromPyObject` over `os.PathLike` and therefore
accepts `pathlib.Path` directly, matching pysam ergonomics.
