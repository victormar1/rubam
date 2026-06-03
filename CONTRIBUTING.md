# Contributing to rubam

Thanks for considering a contribution. This file collects the
non-obvious things you need to know before writing code for rubam.

## Quick start

```bash
# Clone + build (Windows, Linux, macOS — same recipe)
git clone https://github.com/victormar1/rubam
cd rubam

# Use a virtualenv: `maturin develop` requires one.
python -m venv .venv
# Linux/macOS:  source .venv/bin/activate
# Windows:      .venv\Scripts\activate
pip install maturin pytest numpy

maturin develop --release      # builds + installs the extension into the venv

# Test
python -m pytest                          # Python test suite
cargo test --release --no-default-features  # pure-Rust `api` tests
```

`numpy` is only needed for the `get_depths_numpy` return path; everything
else works without it.

## Project structure

```
src/
├── lib.rs                       pyo3 module entry
├── alignment.rs                 BAM read/write pyclasses
├── api/                         Public pure-Rust crate API (no pyo3)
│   ├── mod.rs
│   ├── alignment_file.rs
│   ├── aligned_segment.rs
│   ├── aux_data.rs              ★ NOT aux.rs — see below
│   ├── cigar.rs
│   ├── error.rs
│   └── header.rs
├── common.rs                    shared BAM I/O + tolerant header reader
├── variant.rs                   VCF/BCF read+write pyclasses
├── tools/                       Rust ports of samtools/bcftools subcommands
│   ├── samtools side: sort, index, view, merge, flagstat, idxstats,
│   │   calmd, faidx
│   └── bcftools/: view, norm, concat, query, index, sort, stats
└── bin/                         Standalone shadow CLI binaries
    ├── samtools.rs
    └── bcftools.rs

tests/                           Rust integration tests + Python tests
fuzz/                            cargo-fuzz targets
docs/                            Compatibility / conformance matrices
```

## Windows-reserved filenames — read this first

The following filenames cannot be used on Windows because they are
reserved at the kernel level by Win32:

```
aux  con  prn  nul
com1 com2 com3 com4 com5 com6 com7 com8 com9
lpt1 lpt2 lpt3 lpt4 lpt5 lpt6 lpt7 lpt8 lpt9
```

This applies even with extensions: `aux.rs` cannot be `git add`-ed on
Windows MSVC. **Workaround:** rename the file with a suffix or alternate
name; the exported Rust type can keep its original name. We did this for
`aux.rs` → `aux_data.rs`; the exported type is still `pub struct Aux`.
If you add a file with one of the reserved stems, CI fails on the Windows
runner.

## Standalone CLI binaries

The shadow CLIs (`rubam-samtools`, `rubam-bcftools`, `rubam-depth`,
`rubam-synth-bam`) must be built **without** the default `python` feature
so they do not link pyo3 / the Python runtime (no `python3.dll` dependency
on Windows):

```bash
cargo build --release --bins --no-default-features
```

## Build gotcha (Conda + venv)

If both a Conda environment and a project venv are active, maturin errors
with `Both VIRTUAL_ENV and CONDA_PREFIX are set.` Unset Conda first:

```bash
unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL CONDA_PYTHON_EXE
maturin develop --release
```

## Test fixtures

- `tests/example.bam` — small BAM used for integration tests, committed.
- `tests/fixtures/smoke.bam` — bundled BAM exercised by the write-path and
  pysam-parity tests, committed.
- `tests/data/harmos_compat.bam` — fixture for the public Rust API contract
  tests, committed.
- `tests/data/validation_3sample_100rec.vcf.gz` — synthetic 3-sample VCF
  for the cross-tool validation tests, committed.

Optional large datasets (e.g. NA12878 CRAM + GRCh38 reference) are not
bundled; point `RUBAM_TEST_CRAM` / `RUBAM_TEST_REF` at local copies to
exercise the CRAM path, otherwise those tests skip.

## Pull requests

- Branch off `main`.
- Run `cargo test --no-default-features`, `cargo fmt`, and `pytest` before
  submitting.
- Add a CHANGELOG.md entry (Unreleased section) for any user-visible change.
- For breaking changes, document the migration path in CHANGELOG.md.

## Code style

- Rust: `cargo fmt` defaults.
- Python: 4-space indent, type hints encouraged.
- Commit messages: imperative, ≤ 72 chars first line. Body explains the
  why, not the what.

## Versioning

semver from v0.3.0 onward. v0.x can ship breaking changes between minor
versions but each must be flagged in CHANGELOG.md with a migration path.

## License

MIT. By contributing you agree to license your contribution under MIT.
