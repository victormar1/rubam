# Wheel smoke test

**Last updated: 2026-05-14 (against rubam v0.3.2).**

## Purpose

The `wheel-smoke-test.yml` workflow proves that the **built wheel** is
functional in a clean environment, not just the source tree on a
developer's machine.

The standard CI (`integration.yaml`) runs `maturin develop` and then
`pytest` *inside the source checkout*. That validates the code, but it
does not catch packaging regressions where the wheel is missing a file,
ships an unusable extension module, or fails to import in a fresh venv.

This workflow adds that missing layer:

1. **`build` job** — uses `PyO3/maturin-action` to build the wheel on
   each target OS. Since `rubam` uses `pyo3` with `abi3-py38`, *one
   wheel per OS* covers all of Python 3.8 - 3.13.
2. **`smoke` job** (depends on `build`) — on each OS:
   * downloads the matching wheel artifact,
   * creates a brand-new virtualenv via `python -m venv`,
   * installs the wheel with `pip install <wheel>`,
   * runs `python -m rubam.smoke_test` from a scratch directory with
     `PYTHONPATH=""`, so the rubam source tree (if any) cannot shadow
     the installed package.

## What the smoke test does

`rubam/smoke_test.py`:

1. `import rubam` and print `rubam.__version__`.
2. Locate the bundled `tests/fixtures/smoke.bam` via
   `importlib.resources` (with a filesystem fallback).
3. Call `rubam.get_depths(bam_path, "chr1", 100, 500)`.
4. Assert the returned `positions` and `depths` arrays are non-empty
   and that at least one position has non-zero depth.
5. Exit 0 on success.

## OS matrix

| OS               | Wheel target                               |
|------------------|--------------------------------------------|
| `ubuntu-latest`  | `manylinux_x_y_x86_64` (built in a container) |
| `windows-latest` | `win_amd64`                                |
| `macos-latest`   | `macosx_*_arm64` (Apple Silicon)            |
| `macos-13`       | `macosx_*_x86_64` (Intel)                  |

`macos-13` is kept explicitly because `macos-latest` is now ARM-only.
Both macOS wheels are needed if we want to publish Intel + Apple Silicon
binaries.

## Running locally

```bash
# 1. Build a wheel.
maturin build --release --out dist

# 2. Create a clean venv anywhere outside the source tree.
cd /tmp
python -m venv smoke-venv
source smoke-venv/bin/activate   # Windows: .\smoke-venv\Scripts\activate

# 3. Install the wheel and run the smoke test.
pip install /path/to/rubam/dist/rubam-*.whl
cd $TMPDIR   # or any directory that is NOT the rubam source root
PYTHONPATH="" python -m rubam.smoke_test
```

Expected output:

```
rubam version: 0.3.2
smoke bam: /.../tests/fixtures/smoke.bam
sampled 401 positions, max depth = 10
smoke test OK
```

## Fixture

`tests/fixtures/smoke.bam` is a tiny synthetic BAM (50 reads, 5x
coverage over 1 kb on `chr1`, ~2 KB on disk). It is generated with the
`rubam-synth-bam` binary:

```bash
cargo build --release --bin rubam-synth-bam
./target/release/rubam-synth-bam \
    --output tests/fixtures/smoke.bam \
    --chrom chr1 --length 1000 \
    --coverage 5 --read-length 100 --seed 42
```

Both `smoke.bam` and the auto-generated `smoke.bam.bai` index are
committed so the workflow does not need to regenerate them.

## CI status

The workflow YAML at `.github/workflows/wheel-smoke-test.yml` runs the
4-OS matrix on every push to `main` and on every release tag. It has
**not yet been exercised on GitHub Actions** because the v0.3.x
development happens on a fork that is not pushed to a GitHub remote
(local-only tags per the project's release policy). What has been
verified:

- **Locally on Windows** (the only OS the author has access to): the
  full smoke pipeline runs end-to-end — `maturin develop --release` →
  fresh `python -m venv` outside the source tree → `pip install` of
  the just-built wheel → `python -m rubam.smoke_test` reports
  `rubam version: 0.3.2 / sampled 401 positions, max depth = 10 / smoke test OK`.
- **YAML syntax + step graph** validated by parsing the workflow file
  (`python -c "import yaml; yaml.safe_load(...)"`).
- The job graph (build → smoke with `needs:`) is reviewed for
  correctness; the smoke step explicitly `cd`s out of the source tree
  before running, so the installed wheel is the only `rubam` on the
  Python import path.

The macOS-arm64 and macOS-13-x86_64 cells of the matrix have **not**
been verified locally either — the author has no Mac. They will be
verified the first time this repository is pushed to GitHub and the
workflow runs.

## Windows linkage status of the CLI wrapper binaries

The `rubam-samtools`, `rubam-bcftools`, `rubam-depth` and
`rubam-synth-bam` binaries are built with
`cargo build --release --no-default-features`, which excludes
the optional `python` feature (and therefore the pyo3 / numpy
dependencies) from the binary image. They run on stock Windows
without `python3.dll` on `PATH`. The Python extension module
(`rubam._rubam`) is built with the default features ON via
`maturin build` / `maturin develop`. Pytest suite on Windows 11
Pro: 189 passed, 1 skipped, 2 xfailed (the xfails are the
documented unsupported-codec CRAM cases).

## `pyproject.toml` wheel-include (applied in v0.3.2)

For the smoke test to find `smoke.bam` from an installed wheel, the
fixture must be shipped inside the wheel. The required `include`
directive under `[tool.maturin]` was applied in v0.3.2:

```toml
[tool.maturin]
bindings = "pyo3"
module-name = "rubam._rubam"
features = ["pyo3/extension-module", "pyo3/abi3-py38"]
include = [
    { path = "tests/fixtures/smoke.bam",     format = "wheel" },
    { path = "tests/fixtures/smoke.bam.bai", format = "wheel" },
]
```

Verified locally: `maturin develop --release` produces a wheel
carrying `tests/fixtures/smoke.bam` and `.bai`, and
`python -m rubam.smoke_test` runs against the installed package
reporting `rubam version: 0.3.2 ... smoke test OK`.
