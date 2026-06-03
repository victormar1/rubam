import os
import subprocess
import sys
from pathlib import Path


def _find_rubam_samtools_binary() -> str:
    """Locate the rubam-samtools binary across platforms.

    Order:
      1. RUBAM_SAMTOOLS_BIN env var (for CI / unusual setups).
      2. Windows: target/release/rubam-samtools.exe relative to the repo.
      3. Linux (incl. WSL): ~/.rubam-target-linux/release/rubam-samtools
         (the path setup_wsl.sh writes to).
      4. ./target/release/rubam-samtools (any host).
    """
    env = os.environ.get("RUBAM_SAMTOOLS_BIN")
    if env:
        return env
    repo = Path(__file__).resolve().parent.parent
    suffix = ".exe" if sys.platform == "win32" else ""
    candidates = [
        repo / "target" / "release" / f"rubam-samtools{suffix}",
        Path.home() / ".rubam-target-linux" / "release" / f"rubam-samtools{suffix}",
    ]
    for c in candidates:
        if c.exists():
            return str(c)
    # Last resort: assume PATH has it.
    return f"rubam-samtools{suffix}"


EXE = _find_rubam_samtools_binary()
EXAMPLE_BAM = str(Path(__file__).parent / "example.bam")


def run(*args, input=None):
    return subprocess.run([EXE, *args], input=input, capture_output=True, text=True)


def test_samtools_dispatcher_help():
    r = run("--help")
    assert r.returncode == 0, (r.stdout, r.stderr)
    out = (r.stdout + r.stderr).lower()
    assert "subcommand" in out or "usage" in out


def test_samtools_dispatcher_unknown_subcommand_exits_nonzero():
    r = run("not-a-real-cmd")
    assert r.returncode != 0


import os, shutil, tempfile

def test_samtools_sort_then_index():
    with tempfile.TemporaryDirectory() as td:
        bam_in = os.path.join(td, "in.bam")
        shutil.copy(EXAMPLE_BAM, bam_in)
        bam_out = os.path.join(td, "sorted.bam")
        r = run("sort", "-o", bam_out, bam_in)
        assert r.returncode == 0, (r.stdout, r.stderr)
        assert os.path.exists(bam_out) and os.path.getsize(bam_out) > 0
        r = run("index", bam_out)
        assert r.returncode == 0, (r.stdout, r.stderr)
        assert os.path.exists(bam_out + ".bai")


def test_samtools_sort_missing_output_arg_fails():
    r = run("sort", EXAMPLE_BAM)  # no -o
    assert r.returncode != 0


def test_samtools_index_missing_input_arg_fails():
    r = run("index")
    assert r.returncode != 0


def test_samtools_view_count_only():
    r = run("view", "-c", EXAMPLE_BAM)
    assert r.returncode == 0, (r.stdout, r.stderr)
    assert r.stdout.strip().isdigit()
    assert int(r.stdout.strip()) > 0


def test_samtools_view_with_region(tmp_path):
    out = str(tmp_path / "subset.bam")
    r = run("view", "-b", "-o", out, EXAMPLE_BAM, "chr1:999990-1000010")
    assert r.returncode == 0, (r.stdout, r.stderr)
    assert os.path.exists(out) and os.path.getsize(out) > 0


def test_samtools_merge_two_inputs(tmp_path):
    a = str(tmp_path / "a.bam")
    b = str(tmp_path / "b.bam")
    shutil.copy(EXAMPLE_BAM, a)
    shutil.copy(EXAMPLE_BAM, b)
    out = str(tmp_path / "m.bam")
    r = run("merge", out, a, b)
    assert r.returncode == 0, (r.stdout, r.stderr)
    assert os.path.exists(out) and os.path.getsize(out) > 0


def test_samtools_flagstat_outputs_lines():
    r = run("flagstat", EXAMPLE_BAM)
    assert r.returncode == 0, (r.stdout, r.stderr)
    assert "in total" in r.stdout
    # samtools flagstat output has 17+ lines; we just check the body is multi-line
    assert r.stdout.count("\n") >= 10


def test_samtools_idxstats_one_line_per_chrom():
    r = run("idxstats", EXAMPLE_BAM)
    assert r.returncode == 0, (r.stdout, r.stderr)
    lines = [l for l in r.stdout.splitlines() if l.strip()]
    # samtools idxstats: chrom\tlength\tmapped\tunmapped
    assert any(l.startswith("chr1\t") for l in lines)
    # Each line should have 4 tab-separated columns
    for l in lines:
        parts = l.split("\t")
        assert len(parts) == 4, l
        # length, mapped, unmapped should all be integers
        int(parts[1]); int(parts[2]); int(parts[3])


def test_samtools_faidx_subseq(tmp_path):
    fa = tmp_path / "tiny.fa"
    fa.write_text(">chr1\nACGTACGT\n>chr2\nTTTTAAAA\n")
    r = run("faidx", str(fa), "chr1:1-4")
    assert r.returncode == 0, (r.stdout, r.stderr)
    # samtools faidx prints ">chr1:1-4" then the sequence on the next line.
    assert ">chr1:1-4" in r.stdout
    assert "ACGT" in r.stdout


def test_samtools_faidx_index_only_no_region(tmp_path):
    fa = tmp_path / "tiny.fa"
    fa.write_text(">chr1\nACGTACGT\n")
    fai = str(fa) + ".fai"
    if os.path.exists(fai):
        os.remove(fai)
    r = run("faidx", str(fa))
    assert r.returncode == 0, (r.stdout, r.stderr)
    assert os.path.exists(fai) and os.path.getsize(fai) > 0


def test_samtools_calmd_emits_bam_to_stdout(tmp_path):
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\n" + "N" * 250_000_000 + "\n")
    # Build the .fai first
    r = run("faidx", str(fa))
    assert r.returncode == 0, r.stderr
    # calmd emits the BAM to stdout with -b
    r = subprocess.run(
        [EXE, "calmd", "-b", EXAMPLE_BAM, str(fa)],
        capture_output=True,
    )
    assert r.returncode == 0, r.stderr.decode("utf-8", "replace")
    # Output should be a non-empty BAM (BGZF-magic-prefixed)
    assert len(r.stdout) > 0
    assert r.stdout[:4] == b"\x1f\x8b\x08\x04"  # BGZF / gzip magic
