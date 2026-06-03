"""Acceptance tests for the v0.3.3 AlignedSegment write-side setters.

Covers the property setters (query_name, flag, reference_id,
reference_start, mapping_quality, template_length, query_sequence,
query_qualities, cigarstring, cigartuples), the tag setters
(set_tag / remove_tag), the flag-bit setters (set_is_*), and the
synthetic constructor `AlignedSegment(header=...)`.

Each setter test goes through a write -> reopen -> read cycle so that
both the in-memory mutation AND the BAM serialiser are exercised. A
failure here means either the setter mangled the field, or the writer
dropped it on serialisation.

Skipped wholesale when `tests/fixtures/smoke.bam` is missing.
"""
from __future__ import annotations

import os
import tempfile

import pytest

import rubam

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "smoke.bam")


def _has_smoke_bam() -> bool:
    return os.path.exists(SRC)


pytestmark = pytest.mark.skipif(
    not _has_smoke_bam(), reason="smoke.bam fixture missing"
)


def _write_one(tmp_path, mutate) -> str:
    """Read the first record from smoke.bam, mutate it via `mutate(seg)`,
    write it back, return the output path."""
    out = os.path.join(tmp_path, "out.bam")
    bam_in = rubam.AlignmentFile(SRC, "rb")
    seg = next(iter(bam_in))
    mutate(seg)
    bam_out = rubam.AlignmentFile(out, "wb", template=bam_in)
    bam_out.write(seg)
    bam_out.close()
    bam_in.close()
    return out


def _first_record(path: str):
    bam = rubam.AlignmentFile(path, "rb")
    rec = next(iter(bam))
    bam.close()
    return rec


def test_set_query_name_roundtrip(tmp_path):
    out = _write_one(str(tmp_path), lambda s: setattr(s, "query_name", "renamed_read"))
    rec = _first_record(out)
    assert rec.query_name == "renamed_read"


def test_set_flag_and_mapping_quality_roundtrip(tmp_path):
    def mutate(s):
        s.flag = 0x40  # READ1
        s.mapping_quality = 17

    out = _write_one(str(tmp_path), mutate)
    rec = _first_record(out)
    assert rec.flag == 0x40
    assert rec.is_read1 is True
    assert rec.mapping_quality == 17


def test_set_reference_start_zero_based(tmp_path):
    """0-based external <-> 1-based noodles conversion must round-trip."""
    out = _write_one(str(tmp_path), lambda s: setattr(s, "reference_start", 999))
    rec = _first_record(out)
    assert rec.reference_start == 999


def test_set_template_length_roundtrip(tmp_path):
    out = _write_one(str(tmp_path), lambda s: setattr(s, "template_length", -250))
    rec = _first_record(out)
    assert rec.template_length == -250


def test_set_cigarstring_roundtrip(tmp_path):
    """Same length as the source read (smoke.bam reads are 100bp)."""
    bam_in = rubam.AlignmentFile(SRC, "rb")
    first = next(iter(bam_in))
    n = first.query_length
    bam_in.close()
    cig = f"{n}M"
    out = _write_one(str(tmp_path), lambda s: setattr(s, "cigarstring", cig))
    rec = _first_record(out)
    assert rec.cigarstring == cig


def test_set_cigartuples_roundtrip(tmp_path):
    bam_in = rubam.AlignmentFile(SRC, "rb")
    first = next(iter(bam_in))
    n = first.query_length
    bam_in.close()
    # 10S (n-10)M: soft-clip 10bp then a single match block.
    tuples = [(4, 10), (0, n - 10)]
    out = _write_one(str(tmp_path), lambda s: setattr(s, "cigartuples", tuples))
    rec = _first_record(out)
    # cigartuples returns (op, len) ints. Compare by value.
    assert rec.cigartuples == [(4, 10), (0, n - 10)]
    assert rec.cigarstring == f"10S{n - 10}M"


def test_set_query_sequence_and_qualities_roundtrip(tmp_path):
    bam_in = rubam.AlignmentFile(SRC, "rb")
    first = next(iter(bam_in))
    n = first.query_length
    bam_in.close()
    new_seq = "A" * n
    new_quals = [30] * n

    def mutate(s):
        s.query_sequence = new_seq
        s.query_qualities = new_quals

    out = _write_one(str(tmp_path), mutate)
    rec = _first_record(out)
    assert rec.query_sequence == new_seq
    assert list(rec.query_qualities) == new_quals


def test_set_tag_str_int_float_roundtrip(tmp_path):
    def mutate(s):
        s.set_tag("RG", "rg_test")
        s.set_tag("NH", 3)
        s.set_tag("XF", 1.5)

    out = _write_one(str(tmp_path), mutate)
    rec = _first_record(out)
    assert rec.get_tag("RG") == "rg_test"
    assert rec.get_tag("NH") == 3
    # f32 round-trip — exact for 1.5 (representable).
    assert float(rec.get_tag("XF")) == pytest.approx(1.5)


def test_remove_tag_roundtrip(tmp_path):
    def mutate(s):
        s.set_tag("ZZ", "to_be_removed")
        s.set_tag("NH", 5)
        s.remove_tag("ZZ")

    out = _write_one(str(tmp_path), mutate)
    rec = _first_record(out)
    assert rec.has_tag("NH")
    assert not rec.has_tag("ZZ")


def test_set_tag_validates_name_length():
    bam = rubam.AlignmentFile(SRC, "rb")
    seg = next(iter(bam))
    with pytest.raises(ValueError, match="2 chars"):
        seg.set_tag("A", 1)
    with pytest.raises(ValueError, match="2 chars"):
        seg.set_tag("ABC", 1)
    bam.close()


def test_flag_bit_setters_toggle(tmp_path):
    """Each set_is_* helper must toggle exactly its bit (no neighbours)."""
    def mutate(s):
        # Start from a known state — clear everything then toggle in two bits.
        s.flag = 0
        s.set_is_paired(True)
        s.set_is_proper_pair(True)
        s.set_is_read1(True)
        s.set_is_reverse(False)
        s.set_is_duplicate(True)
        # Now turn proper_pair back off — the others must survive.
        s.set_is_proper_pair(False)

    out = _write_one(str(tmp_path), mutate)
    rec = _first_record(out)
    assert rec.is_paired is True
    assert rec.is_proper_pair is False
    assert rec.is_read1 is True
    assert rec.is_reverse is False
    assert rec.is_duplicate is True


def test_constructor_synthesises_record(tmp_path):
    """Pattern 2 from the spec: build a fresh record from scratch."""
    bam_in = rubam.AlignmentFile(SRC, "rb")
    out = os.path.join(tmp_path, "synth.bam")
    bam_out = rubam.AlignmentFile(out, "wb", template=bam_in)

    seg = rubam.AlignedSegment(header=bam_in.header)
    seg.query_name = "synth_1"
    seg.flag = 0
    seg.reference_id = 0
    seg.reference_start = 100
    seg.mapping_quality = 60
    seg.cigarstring = "150M"
    seg.query_sequence = "A" * 150
    seg.query_qualities = [30] * 150
    bam_out.write(seg)
    bam_out.close()
    bam_in.close()

    rec = _first_record(out)
    assert rec.query_name == "synth_1"
    assert rec.reference_start == 100
    assert rec.mapping_quality == 60
    assert rec.cigarstring == "150M"
    assert rec.query_sequence == "A" * 150
    assert list(rec.query_qualities) == [30] * 150


def test_reopen_count_after_synth_write(tmp_path):
    """Verify rubam can reopen its own write and count exactly one record."""
    bam_in = rubam.AlignmentFile(SRC, "rb")
    out = os.path.join(tmp_path, "one.bam")
    bam_out = rubam.AlignmentFile(out, "wb", template=bam_in)

    seg = rubam.AlignedSegment(header=bam_in.header)
    seg.query_name = "only_one"
    seg.flag = 4  # unmapped — avoids index/contig consistency checks
    seg.cigarstring = "100M"
    seg.query_sequence = "C" * 100
    seg.query_qualities = [20] * 100
    bam_out.write(seg)
    bam_out.close()
    bam_in.close()

    bam_re = rubam.AlignmentFile(out, "rb")
    n = sum(1 for _ in iter(bam_re))
    bam_re.close()
    assert n == 1
