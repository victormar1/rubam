# rubam — fuzzing harness

This directory hosts the [`cargo-fuzz`](https://rust-fuzz.github.io/book/cargo-fuzz.html)
infrastructure for `rubam`. The targets exercise the parsers that `rubam`
relies on for its core surfaces:

| Target              | Surface                                                                                  |
| ------------------- | ---------------------------------------------------------------------------------------- |
| `fuzz_bam_reader`   | `noodles::bam::io::Reader` — same code path as `rubam::api::AlignmentFile::open`.        |
| `fuzz_vcf_parser`   | `noodles::vcf::io::Reader` — same code path as `rubam::variant::VariantFile::open`.      |
| `fuzz_cigar_walker` | `noodles::sam::record::Cigar` text parser + an op-walk computing reference/query spans. |
| `fuzz_aux_tags`     | `noodles::sam` aux-tag decoder, the parser feeding `rubam::api::aux_data::aux_from_noodles`. |

Each target follows the libFuzzer contract: arbitrary input bytes go in,
the parser must return `Err` (not panic) on malformed input. A panic, an
abort, or a sanitizer hit counts as a finding.

## Why this is a separate crate

The fuzz crate is not part of the main `rubam` workspace on purpose.

1. `rubam` declares `pyo3 = { features = ["extension-module"] }`.
   On Linux, that feature suppresses the libpython link line, which
   makes the rubam `rlib` unsuitable as a transitive dependency of a
   plain Rust binary like a libFuzzer harness — the link step would
   complain about missing `Py*` symbols.
2. The set of parsers exercised below is exactly the set rubam wraps.
   Fuzzing the `noodles` entry points used by `rubam` is therefore an
   accurate test of rubam's parser robustness without dragging the
   pyo3 build chain into the fuzz harness.

If we ever split rubam into `rubam-core` (pure-Rust) + `rubam-py`
(pyo3 layer), the fuzz crate could depend directly on `rubam-core`.

## Why Linux-only in CI

`cargo-fuzz` is built on libFuzzer's compiler-rt instrumentation, which
is well supported on Linux (`x86_64-unknown-linux-gnu`) and effectively
unsupported on Windows (`x86_64-pc-windows-msvc`) and macOS (the SIP /
codesign story is fragile, and Apple has shipped multiple breaks of
sanitizer support over the past few releases). Upstream
[rust-fuzz/cargo-fuzz#147](https://github.com/rust-fuzz/cargo-fuzz/issues/147)
tracks the Windows side; the consensus is "run fuzz on Linux".

This matches the convention used by every rust-fuzz–style project we
checked (e.g. `serde-rs`, `image-rs`, `regex`). The Windows-native
build / test path remains covered by the main CI workflow.

## Local run (Linux)

```bash
# One-time install. Skipped by default in this repo because the user
# policy bans `cargo install` without explicit ask.
cargo install cargo-fuzz --locked

# List the targets.
cargo fuzz list

# Run one target for 60s.
cargo fuzz run fuzz_bam_reader -- -max_total_time=60

# Reproduce a previously found crash.
cargo fuzz run fuzz_bam_reader fuzz/artifacts/fuzz_bam_reader/crash-<hash>
```

`cargo fuzz` requires the **nightly** toolchain:

```bash
rustup toolchain install nightly
rustup component add --toolchain nightly rust-src
```

## Local run (Windows / macOS)

Not supported by `cargo-fuzz` upstream. Two workarounds for local
investigation:

- **WSL** on Windows: install Ubuntu, then run the Linux command
  above inside WSL. The harness operates on file bytes only — no
  Windows-specific code paths are tested by fuzzing, so WSL is fine.
- **Docker** on macOS: `docker run --rm -it -v "$PWD":/work -w /work
  rust:latest bash -lc 'rustup default nightly && cargo install
  cargo-fuzz --locked && cargo fuzz run fuzz_bam_reader -- -max_total_time=60'`.

## Continuous fuzzing

The nightly schedule lives in
`.github/workflows/fuzz-nightly.yml`. Each target runs for 300 s on a
fresh `ubuntu-latest` runner. Corpora are cached between runs so
coverage accumulates over time. Crashes upload as workflow artifacts.

## Adding a new target

1. Create `fuzz/fuzz_targets/fuzz_<name>.rs` using the
   `#![no_main]` + `libfuzzer_sys::fuzz_target!` skeleton.
2. Add the `[[bin]]` entry to `fuzz/Cargo.toml`.
3. Add the target name to the matrix in
   `.github/workflows/fuzz-nightly.yml`.
