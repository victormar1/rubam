"""Builder/mutation-API acceptance tests for the BAM write path.

The v0.3.2 write path only supports unmodified pass-through records. This
module exercises the v0.3.3 mutation surface added to ``AlignedSegment``:

  1. Mutate an existing record obtained from ``AlignmentFile.fetch(...)``
     (set_flag, set_cigarstring, set_tag), write it back and assert the
     mutations stuck after re-reading.
  2. Synthesise a brand-new ``AlignedSegment`` bound to a Header, set
     every required field, write it to a fresh BAM and roundtrip-verify
     each field.
  3. Cross-validate the synthesised BAM with system samtools through
     WSL when available.

Skipped wholesale if the bundled ``tests/fixtures/smoke.bam`` is missing.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile

import pytest

import rubam
from _wsl_probe import wsl_usable


SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "smoke.bam")


def _has_smoke_bam() -> bool:
    return os.path.exists(SRC)


@pytest.fixture
def tmp_bam(tmp_path):
    return str(tmp_path / "out.bam")


# ---------------------------------------------------------------------------
# (1) Mutate fetched records.
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _has_smoke_bam(), reason="smoke.bam fixture missing")
def test_mutate_fetched_record_flag_cigar_tag_roundtrip(tmp_bam):
    bam_in = rubam.AlignmentFile(SRC, "rb")
    bam_out = rubam.AlignmentFile(tmp_bam, "wb", template=bam_in)

    target_flag = 0x4 | 0x100  # unmapped + secondary (cheap, unambiguous bits)
    # smoke.bam reads are 100bp long; keep the CIGAR consistent so noodles
    # accepts the mutated record.
    target_cigar = "100M"
    target_tag_name = "ZZ"
    target_tag_value = 42

    for r in iter(bam_in):
        r.flag = target_flag
        r.cigarstring = target_cigar
        r.set_tag(target_tag_name, target_tag_value)
        bam_out.write(r)
    bam_out.close()
    bam_in.close()

    bam_re = rubam.AlignmentFile(tmp_bam, "rb")
    seen = 0
    for r in iter(bam_re):
        assert r.flag == target_flag, f"flag not preserved: {r.flag:#x}"
        assert r.cigarstring == target_cigar, f"cigar: {r.cigarstring!r}"
        assert r.has_tag(target_tag_name), "ZZ tag missing after roundtrip"
        assert r.get_tag(target_tag_name) == target_tag_value
        seen += 1
    bam_re.close()
    assert seen == 50, f"expected 50 records, got {seen}"


@pytest.mark.skipif(not _has_smoke_bam(), reason="smoke.bam fixture missing")
def test_mutate_query_name_and_sequence_roundtrip(tmp_bam):
    bam_in = rubam.AlignmentFile(SRC, "rb")
    bam_out = rubam.AlignmentFile(tmp_bam, "wb", template=bam_in)

    new_seq = "ACGTACGTAC"
    new_quals = [30, 31, 32, 33, 34, 35, 36, 37, 38, 39]

    first = True
    expected_first_name = None
    for r in iter(bam_in):
        if first:
            expected_first_name = "RENAMED_READ_0"
            r.query_name = expected_first_name
            r.query_sequence = new_seq
            r.query_qualities = new_quals
            # Make cigar consistent with the new 10-base sequence so noodles
            # can encode it without complaint.
            r.cigarstring = "10M"
            first = False
        bam_out.write(r)
    bam_out.close()
    bam_in.close()

    bam_re = rubam.AlignmentFile(tmp_bam, "rb")
    first_after = next(iter(bam_re))
    assert first_after.query_name == expected_first_name
    assert first_after.query_sequence == new_seq
    assert list(first_after.query_qualities) == new_quals
    assert first_after.cigarstring == "10M"
    bam_re.close()


# ---------------------------------------------------------------------------
# (2) Synthesise from scratch.
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _has_smoke_bam(), reason="smoke.bam fixture missing")
def test_synthesise_fresh_aligned_segment_and_write(tmp_bam):
    src = rubam.AlignmentFile(SRC, "rb")
    hdr = src.header
    src.close()

    bam_out = rubam.AlignmentFile(tmp_bam, "wb", header=hdr)
    seg = rubam.AlignedSegment(hdr)
    seg.query_name = "synth_read_1"
    seg.flag = 0  # mapped, primary, single-end
    seg.reference_id = 0
    seg.reference_start = 99  # 0-based -> pysam-like
    seg.mapping_quality = 60
    seg.cigarstring = "5M"
    seg.query_sequence = "ACGTA"
    seg.query_qualities = [40, 41, 42, 43, 44]
    seg.template_length = 0
    seg.set_tag("NM", 0)
    seg.set_tag("XS", "rubam-synth")

    bam_out.write(seg)
    bam_out.close()

    bam_re = rubam.AlignmentFile(tmp_bam, "rb")
    recs = list(iter(bam_re))
    bam_re.close()

    assert len(recs) == 1
    got = recs[0]
    assert got.query_name == "synth_read_1"
    assert got.flag == 0
    assert got.reference_id == 0
    assert got.reference_start == 99
    assert got.mapping_quality == 60
    assert got.cigarstring == "5M"
    assert got.query_sequence == "ACGTA"
    assert list(got.query_qualities) == [40, 41, 42, 43, 44]
    assert got.template_length == 0
    assert got.has_tag("NM") and got.get_tag("NM") == 0
    assert got.has_tag("XS") and got.get_tag("XS") == "rubam-synth"


@pytest.mark.skipif(not _has_smoke_bam(), reason="smoke.bam fixture missing")
def test_synthesise_via_cigartuples(tmp_bam):
    """The cigartuples setter parallel to cigarstring."""
    src = rubam.AlignmentFile(SRC, "rb")
    hdr = src.header
    src.close()

    bam_out = rubam.AlignmentFile(tmp_bam, "wb", header=hdr)
    seg = rubam.AlignedSegment(hdr)
    seg.query_name = "ct_synth"
    seg.flag = 0
    seg.reference_id = 0
    seg.reference_start = 0
    seg.mapping_quality = 30
    # 3M1I2M, i.e. matches + 1bp insertion + 2 more matches; read length 6.
    seg.cigartuples = [(0, 3), (1, 1), (0, 2)]
    seg.query_sequence = "ACGTAC"
    seg.query_qualities = [20, 21, 22, 23, 24, 25]
    bam_out.write(seg)
    bam_out.close()

    bam_re = rubam.AlignmentFile(tmp_bam, "rb")
    got = next(iter(bam_re))
    bam_re.close()
    assert got.cigartuples == [(0, 3), (1, 1), (0, 2)]
    assert got.cigarstring == "3M1I2M"


# ---------------------------------------------------------------------------
# (3) set_tags bulk + to_dict / from_dict roundtrip.
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _has_smoke_bam(), reason="smoke.bam fixture missing")
def test_set_tags_bulk_replaces_existing(tmp_bam):
    src = rubam.AlignmentFile(SRC, "rb")
    hdr = src.header
    src.close()

    bam_out = rubam.AlignmentFile(tmp_bam, "wb", header=hdr)
    seg = rubam.AlignedSegment(hdr)
    seg.query_name = "tagged"
    seg.flag = 4  # unmapped — no cigar/pos requirements
    seg.set_tag("OL", 999)  # will be wiped by the bulk set_tags below
    seg.tags = [("AA", 1), ("BB", "two"), ("CC", 3.5)]
    bam_out.write(seg)
    bam_out.close()

    bam_re = rubam.AlignmentFile(tmp_bam, "rb")
    got = next(iter(bam_re))
    bam_re.close()
    names = {t[0] for t in got.tags}
    assert names == {"AA", "BB", "CC"}, f"set_tags did not replace OL: {names}"
    assert got.get_tag("AA") == 1
    assert got.get_tag("BB") == "two"
    assert abs(got.get_tag("CC") - 3.5) < 1e-6


@pytest.mark.skipif(not _has_smoke_bam(), reason="smoke.bam fixture missing")
def test_to_dict_from_dict_roundtrip(tmp_bam):
    """to_dict() emits the format that from_dict() can ingest."""
    src = rubam.AlignmentFile(SRC, "rb")
    hdr = src.header
    first = next(iter(src))
    d = first.to_dict()
    src.close()

    # Sanity-check the dict surface.
    for key in ("name", "flag", "ref_name", "ref_pos", "map_quality",
                "cigar", "next_ref_name", "next_ref_pos", "length",
                "seq", "qual", "tags"):
        assert key in d, f"to_dict missing key {key!r}"

    # Rebuild from the dict, write, re-read, compare scalar fields.
    seg = rubam.AlignedSegment.from_dict(hdr, d)
    bam_out = rubam.AlignmentFile(tmp_bam, "wb", header=hdr)
    bam_out.write(seg)
    bam_out.close()

    bam_re = rubam.AlignmentFile(tmp_bam, "rb")
    got = next(iter(bam_re))
    bam_re.close()

    assert got.query_name == d["name"]
    assert got.flag == d["flag"]
    # cigar may be None if the source record had no CIGAR.
    assert got.cigarstring == d["cigar"]
    assert got.mapping_quality == d["map_quality"]
    if d["ref_pos"] is not None:
        assert got.reference_start == d["ref_pos"] - 1


# ---------------------------------------------------------------------------
# (4) Cross-validate with system samtools (WSL on Windows).
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    sys.platform != "win32" or not wsl_usable(),
    reason="WSL Ubuntu (system samtools) not usable — cross-validation skipped",
)
def test_synth_record_validates_against_samtools(tmp_bam):
    """End-to-end: build a single synthetic record from scratch, write
    it, and verify samtools accepts the BAM and counts 1 record."""
    if not _has_smoke_bam():
        pytest.skip("smoke.bam missing (needed for the header template)")
    src = rubam.AlignmentFile(SRC, "rb")
    hdr = src.header
    src.close()

    bam_out = rubam.AlignmentFile(tmp_bam, "wb", header=hdr)
    seg = rubam.AlignedSegment(hdr)
    seg.query_name = "samtools_xval"
    seg.flag = 0
    seg.reference_id = 0
    seg.reference_start = 0
    seg.mapping_quality = 30
    seg.cigarstring = "4M"
    seg.query_sequence = "ACGT"
    seg.query_qualities = [30, 31, 32, 33]
    bam_out.write(seg)
    bam_out.close()

    drive = tmp_bam[0].lower()
    wsl_path = "/mnt/" + drive + tmp_bam[2:].replace("\\", "/")
    proc = subprocess.run(
        ["wsl", "-d", "Ubuntu", "samtools", "view", "-c", wsl_path],
        capture_output=True, text=True, timeout=30,
    )
    assert proc.returncode == 0, f"samtools rejected the BAM: {proc.stderr}"
    assert proc.stdout.strip() == "1", (
        f"samtools view -c returned {proc.stdout!r}, expected '1'"
    )
