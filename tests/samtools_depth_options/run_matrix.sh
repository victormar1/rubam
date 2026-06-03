#!/usr/bin/env bash
# tests/samtools_depth_options/run_matrix.sh
#
# Empirical compatibility matrix between `samtools depth` (system) and
# `rubam-samtools depth`. Runs both for each option combo, diffs stdout,
# classifies the result. Output: ./results/ + a printable summary.
#
# Designed to run on WSL Ubuntu (default distro) against the project
# binaries built in target/release/. The Windows .exe files are invoked
# directly thanks to WSL interop.
#
# Usage:
#     bash tests/samtools_depth_options/run_matrix.sh
#
# Exit code: 0 if the matrix ran end-to-end (regardless of how many
# options diverge), 2 if a hard prerequisite is missing (samtools or
# rubam binaries).
set -u
set -o pipefail

# --- locate project root ---------------------------------------------------
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"

RESULTS_DIR="$HERE/results"
mkdir -p "$RESULTS_DIR"

# --- prerequisites ---------------------------------------------------------
if ! command -v samtools >/dev/null 2>&1; then
    echo "FATAL: system samtools not found in PATH" >&2
    exit 2
fi

# Pick the rubam-samtools binary: release first, then debug. Accept either
# the bare unix name (Linux build) or the .exe (Windows build via WSL interop).
pick_bin() {
    local name="$1"
    for cand in \
        "$ROOT/target/release/$name" \
        "$ROOT/target/release/$name.exe" \
        "$ROOT/target/debug/$name" \
        "$ROOT/target/debug/$name.exe"; do
        if [[ -x "$cand" ]]; then
            echo "$cand"
            return 0
        fi
    done
    return 1
}

RUBAM_SAMTOOLS="$(pick_bin rubam-samtools || true)"
RUBAM_DEPTH="$(pick_bin rubam-depth || true)"
if [[ -z "${RUBAM_SAMTOOLS:-}" ]]; then
    echo "FATAL: rubam-samtools binary not found. Run 'maturin develop' or 'cargo build --release' first." >&2
    exit 2
fi
if [[ -z "${RUBAM_DEPTH:-}" ]]; then
    echo "WARN: rubam-depth binary not found — the 'fallback to rubam-depth' classification will be skipped." >&2
fi

# When the rubam binary is a Windows .exe and we are inside WSL, paths must
# be translated to Win32 form (D:\...) before being passed as argv. samtools
# is a native Linux binary and keeps using POSIX paths.
IS_WIN_BIN=0
if [[ "$RUBAM_SAMTOOLS" == *.exe ]] && command -v wslpath >/dev/null 2>&1; then
    IS_WIN_BIN=1
fi
to_native() {
    # Convert a POSIX path to Win32 if the rubam binary is a .exe, else echo as-is.
    if [[ "$IS_WIN_BIN" == 1 ]]; then
        wslpath -w "$1"
    else
        echo "$1"
    fi
}

# --- fixture BAM -----------------------------------------------------------
# Order: a small fixture in tests/, else the 9.8 MB synthetic chr20 BAM.
FIXTURE=""
for cand in \
    "$ROOT/tests/fixtures/depth_fixture.bam" \
    "$ROOT/data/synthetic/synth_chr20_1Mb_30x.bam" \
    "$ROOT/tests/example.bam"; do
    if [[ -f "$cand" ]]; then
        FIXTURE="$cand"
        break
    fi
done
if [[ -z "$FIXTURE" ]]; then
    echo "FATAL: no fixture BAM found" >&2
    exit 2
fi

# Make sure the fixture has an index — samtools depth needs it for -r, and
# rubam-depth uses an indexed reader unconditionally.
if [[ ! -f "${FIXTURE}.bai" && ! -f "${FIXTURE%.bam}.bai" ]]; then
    echo "[setup] indexing $FIXTURE"
    samtools index "$FIXTURE"
fi

# Use the first @SQ name as the canonical chromosome for region tests.
CHROM="$(samtools view -H "$FIXTURE" | awk '$1=="@SQ" {for(i=2;i<=NF;i++) if($i ~ /^SN:/){sub("^SN:","",$i); print $i; exit}}')"
CHROM="${CHROM:-chr20}"

# A short region every test reuses so diffs stay tiny.
REGION="${CHROM}:1-1000"
REGION_START=1
REGION_END=1000

# A tiny BED file: same region.
BED="$RESULTS_DIR/region.bed"
printf '%s\t%d\t%d\n' "$CHROM" 0 "$REGION_END" > "$BED"

# Path variants used when invoking the rubam binaries (Win32 form under WSL).
FIXTURE_NATIVE="$(to_native "$FIXTURE")"
BED_NATIVE="$(to_native "$BED")"
RESULTS_DIR_NATIVE="$(to_native "$RESULTS_DIR")"

echo "[setup] fixture = $FIXTURE"
echo "[setup] chrom   = $CHROM"
echo "[setup] region  = $REGION"
echo "[setup] samtools: $(samtools --version | head -1)"
echo "[setup] rubam-samtools: $RUBAM_SAMTOOLS"
echo "[setup] rubam-depth   : ${RUBAM_DEPTH:-<missing>}"
echo

# --- harness ---------------------------------------------------------------
# A single test row.  Args:
#   $1 = option label (e.g. "-a", "-r REGION")
#   $2 = samtools argv (string, will be eval'd)
#   $3 = rubam-samtools argv suffix (string), or "__SKIP__" to skip
#   $4 = optional alternate rubam-depth argv (string), or empty
#   $5 = notes (free text)
#
# Writes:
#   $RESULTS_DIR/<slug>.samtools.tsv
#   $RESULTS_DIR/<slug>.rubam.tsv         (only if rubam-samtools accepts)
#   $RESULTS_DIR/<slug>.rubamdepth.tsv    (only if rubam-depth was invoked)
#   appends a row to $RESULTS_DIR/summary.tsv
slugify() {
    echo "$1" | tr ' /:-' '____' | tr -cd 'A-Za-z0-9_'
}

SUMMARY="$RESULTS_DIR/summary.tsv"
printf 'option\tsupported\toutput_match\tnotes\n' > "$SUMMARY"

run_one() {
    local label="$1" samtools_args="$2" rubam_args="$3" rubam_depth_args="${4:-}" notes="${5:-}"
    local slug; slug="$(slugify "$label")"
    local s_out="$RESULTS_DIR/${slug}.samtools.tsv"
    local r_out="$RESULTS_DIR/${slug}.rubam.tsv"
    local d_out="$RESULTS_DIR/${slug}.rubamdepth.tsv"
    local s_err="$RESULTS_DIR/${slug}.samtools.err"
    local r_err="$RESULTS_DIR/${slug}.rubam.err"
    local d_err="$RESULTS_DIR/${slug}.rubamdepth.err"

    echo "--- $label ---"

    # 1. system samtools (reference): POSIX paths.
    eval "samtools depth $samtools_args \"\$FIXTURE\"" > "$s_out" 2> "$s_err" || true

    # 2. rubam-samtools depth (the headline claim): Win32 paths under WSL.
    local supported="no" outmatch="n/a"
    if [[ "$rubam_args" != "__SKIP__" ]]; then
        local rubam_args_native="${rubam_args//$BED/$BED_NATIVE}"
        rubam_args_native="${rubam_args_native//$RESULTS_DIR/$RESULTS_DIR_NATIVE}"
        eval "\"\$RUBAM_SAMTOOLS\" depth $rubam_args_native \"\$FIXTURE_NATIVE\"" > "$r_out" 2> "$r_err" || true
        if [[ -s "$r_out" ]] && ! grep -q 'unknown subcommand' "$r_err"; then
            if diff -q "$s_out" "$r_out" >/dev/null 2>&1; then
                supported="yes"; outmatch="exact"
            elif diff -q <(tr -s '[:space:]' ' ' < "$s_out") \
                          <(tr -s '[:space:]' ' ' < "$r_out") >/dev/null 2>&1; then
                supported="yes"; outmatch="byte-equivalent"
            else
                supported="partial"; outmatch="diverges"
            fi
        else
            supported="no"; outmatch="n/a"
        fi
    fi

    # 3. fallback probe: does rubam-depth implement an equivalent?
    local depth_fallback=""
    if [[ -n "$rubam_depth_args" && -n "${RUBAM_DEPTH:-}" ]]; then
        eval "\"\$RUBAM_DEPTH\" \"\$FIXTURE_NATIVE\" $rubam_depth_args" > "$d_out" 2> "$d_err" || true
        if [[ -s "$d_out" ]]; then
            if diff -q "$s_out" "$d_out" >/dev/null 2>&1; then
                depth_fallback="rubam-depth: exact"
            elif diff -q <(tr -s '[:space:]' ' ' < "$s_out") \
                          <(tr -s '[:space:]' ' ' < "$d_out") >/dev/null 2>&1; then
                depth_fallback="rubam-depth: byte-equivalent"
            else
                depth_fallback="rubam-depth: diverges"
            fi
        else
            depth_fallback="rubam-depth: error ($(head -1 "$d_err" 2>/dev/null | tr -d '\r' | head -c 80))"
        fi
    fi

    local full_notes="$notes"
    if [[ -n "$depth_fallback" ]]; then
        if [[ -n "$full_notes" ]]; then
            full_notes="$full_notes; $depth_fallback"
        else
            full_notes="$depth_fallback"
        fi
    fi
    if [[ "$supported" == "no" && -s "$r_err" ]]; then
        local err_msg; err_msg="$(head -1 "$r_err" | tr -d '\r' | head -c 100)"
        if [[ -n "$err_msg" ]]; then
            if [[ -n "$full_notes" ]]; then
                full_notes="rubam-samtools err: \"$err_msg\"; $full_notes"
            else
                full_notes="rubam-samtools err: \"$err_msg\""
            fi
        fi
    fi

    printf '%s\t%s\t%s\t%s\n' "$label" "$supported" "$outmatch" "$full_notes" >> "$SUMMARY"
    echo "  -> supported=$supported  match=$outmatch  $full_notes"
}

# --- matrix ----------------------------------------------------------------
# Baseline (no flag): samtools default = skip zero-depth; rubam-depth always
# emits every position. We test this first because the divergence on the
# default is itself a documented finding.
run_one "(no flag, baseline)" \
        "-r $REGION" \
        "-r $REGION" \
        "$CHROM $REGION_START $REGION_END" \
        "samtools default skips zero-depth positions; rubam-depth emits all"

run_one "-a" \
        "-a -r $REGION" \
        "-a -r $REGION" \
        "$CHROM $REGION_START $REGION_END" \
        "emit zero-depth positions within the region"

run_one "-aa" \
        "-aa -r $REGION" \
        "-aa -r $REGION" \
        "$CHROM $REGION_START $REGION_END" \
        "absolutely all positions, incl. unused refs"

run_one "-q 13" \
        "-a -q 13 -r $REGION" \
        "-a -q 13 -r $REGION" \
        "$CHROM $REGION_START $REGION_END -q 13" \
        "min base quality"

run_one "-Q 10" \
        "-a -Q 10 -r $REGION" \
        "-a -Q 10 -r $REGION" \
        "$CHROM $REGION_START $REGION_END -Q 10" \
        "min mapping quality"

run_one "-r REGION" \
        "-a -r $REGION" \
        "-a -r $REGION" \
        "$CHROM $REGION_START $REGION_END" \
        "region restriction"

run_one "-b BED" \
        "-a -b $BED" \
        "-a -b $BED" \
        "" \
        "BED file region restriction"

run_one "-G 0x4" \
        "-a -G 0x4 -r $REGION" \
        "-a -G 0x4 -r $REGION" \
        "" \
        "additional filter-out flags"

run_one "-d 100" \
        "-a -d 100 -r $REGION" \
        "-a -d 100 -r $REGION" \
        "$CHROM $REGION_START $REGION_END -d 100" \
        "max depth cap"

run_one "-H (header line)" \
        "-a -H -r $REGION" \
        "-a -H -r $REGION" \
        "" \
        "print column header line"

run_one "-o FILE" \
        "-a -o $RESULTS_DIR/_dummy_samtools.tsv -r $REGION; cat $RESULTS_DIR/_dummy_samtools.tsv" \
        "-a -o $RESULTS_DIR/_dummy_rubam.tsv -r $REGION; cat $RESULTS_DIR/_dummy_rubam.tsv 2>/dev/null" \
        "" \
        "write output to FILE"

# --- print summary ---------------------------------------------------------
echo
echo "================= MATRIX SUMMARY ================="
column -t -s $'\t' "$SUMMARY"
echo "=================================================="
echo "Per-test artifacts: $RESULTS_DIR/"
echo "Summary TSV       : $SUMMARY"
