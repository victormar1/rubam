# Security Policy

## Reporting a Vulnerability

If you discover a security vulnerability in **rubam**, please report it
privately by email to:

- **civi64@gmail.com**

Please include:

- A description of the issue and its impact.
- Steps to reproduce (input file, command line, version of rubam, OS).
- Any proof-of-concept code or sample BAM/VCF/BCF/CRAM that triggers the bug.
- Whether you intend to disclose publicly and on what timeline.

Do **not** open a public GitHub issue for security-sensitive reports.

### Response SLA

- We will acknowledge receipt within **14 calendar days**.
- We will provide a triage assessment (accepted / needs-info / not-a-vuln)
  within the same window.
- Fix timelines depend on severity; critical issues are prioritized over
  feature work.

Coordinated disclosure is welcome — propose an embargo date and we will
work with you on a coordinated release.

## Supported Versions

Security fixes are backported to currently supported minor lines only.

| Version | Supported          |
| ------- | ------------------ |
| 0.4.x   | Yes (once released) |
| 0.3.x   | Yes                |
| < 0.3   | No                 |

When a new minor (e.g. 0.5.x) ships, the oldest supported line is dropped
from this table.

## Supply-Chain Hardening

rubam enforces the following on every PR and on a weekly cron:

- **`cargo audit`** — RustSec advisory database. Blocks merges on any known
  vulnerability in the dependency graph.
- **`cargo deny check`** — license allowlist, source allowlist (crates.io
  only), duplicate-major-version warnings, advisory deny. Config lives in
  [`deny.toml`](./deny.toml).
- **CycloneDX SBOM** — generated on each release tag (`v*`) via
  `cargo cyclonedx` and uploaded as a workflow artifact for downstream
  attestation. See `.github/workflows/sbom.yml`.
- **Dependabot** — weekly updates for both the `cargo` and `github-actions`
  ecosystems. Config in `.github/dependabot.yml`.

### Reproducing the checks locally

The supply-chain tooling is **not** vendored. Install on demand:

```bash
# Audit (RustSec advisories)
cargo install cargo-audit --locked
cargo audit

# Deny (licenses, bans, sources, advisories)
cargo install cargo-deny --locked
cargo deny check

# SBOM (CycloneDX JSON + XML)
cargo install cargo-cyclonedx --locked
cargo cyclonedx --format json
cargo cyclonedx --format xml
```

All three are pinned via `--locked` to honor `Cargo.lock` and avoid silent
transitive upgrades during install.

## Threat Model — what is in scope

- Memory-safety bugs reachable from user-supplied BAM/VCF/BCF/CRAM input.
- Path-traversal or symlink-escape in index/output file handling.
- Denial-of-service via pathological compressed inputs (zip-bomb-style).
- Vulnerabilities in pinned dependencies that affect rubam's API surface.

## Out of scope

- Bugs requiring an attacker-controlled `Cargo.toml` (i.e. the user must
  already have full code execution on the build host).
- Issues in `pyo3` extension loading on systems with broken Python ABIs —
  report those upstream to PyO3.
- Performance regressions without a security impact.
