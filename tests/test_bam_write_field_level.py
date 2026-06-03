"""Field-level roundtrip validation for the BAM write path.

The reviewer's v5 critique: `samtools view -c` only proves the file is
readable and has the right record count. It does NOT prove that flags,
CIGAR, tags, qualities, template_length, mate references, etc. survive
the read -> write -> re-read cycle intact.

This module closes that gap. For every record in the source fixture we
read it, write it via `rubam.AlignmentFile.write()`, then re-read the
output file and compare every exposed field on the rubam side. A
failure here means the BAM writer dropped or mangled a field.

Skipped if `tests/fixtures/smoke.bam` is missing (it is bundled with the
wheel from v0.3.2 onwards, so it should be present).
"""
from __future__ import annotations

import os
import tempfile

import pytest

import rubam


SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "smoke.bam")

FIELDS = (
    "query_name",
    "flag",
    "reference_id",
    "reference_name",
    "reference_start",
    "reference_end",
    "mapping_quality",
    "cigarstring",
    "cigartuples",
    "query_sequence",
    "query_qualities",
    "query_length",
    "template_length",
    "is_paired",
    "is_proper_pair",
    "is_unmapped",
    "is_mate_unmapped",
    "is_reverse",
    "is_mate_reverse",
    "is_read1",
    "is_read2",
    "is_secondary",
    "is_qcfail",
    "is_duplicate",
    "is_supplementary",
)


def _snapshot(record) -> dict:
    """Materialise every field of a record into a plain dict for comparison."""
    snap = {}
    for f in FIELDS:
        try:
            v = getattr(record, f)
        except Exception as e:
            v = f"<error: {type(e).__name__}: {e}>"
        snap[f] = v
    try:
        tags = sorted(list(record.tags), key=lambda nv: nv[0])
        snap["_tags"] = tags
    except Exception as e:
        snap["_tags"] = f"<error: {type(e).__name__}>"
    return snap


@pytest.mark.skipif(not os.path.exists(SRC), reason="smoke.bam fixture missing")
def test_field_level_roundtrip_every_record():
    """Read N records, snapshot them, write them, re-read, snapshot again,
    assert every field matches record-by-record."""
    bam_in = rubam.AlignmentFile(SRC, "rb")
    before = [_snapshot(r) for r in iter(bam_in)]
    bam_in.close()
    assert before, "fixture must contain at least one record"

    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "roundtrip.bam")
        bam_in2 = rubam.AlignmentFile(SRC, "rb")
        bam_out = rubam.AlignmentFile(out, "wb", template=bam_in2)
        for r in iter(bam_in2):
            bam_out.write(r)
        bam_out.close()
        bam_in2.close()

        rubam.index(out)
        bam_re = rubam.AlignmentFile(out, "rb")
        after = [_snapshot(r) for r in iter(bam_re)]
        bam_re.close()

    assert len(after) == len(before), (
        f"record count differs: read {len(before)}, wrote-then-read {len(after)}"
    )
    for i, (b, a) in enumerate(zip(before, after)):
        for f in FIELDS:
            assert b[f] == a[f], (
                f"record #{i} field {f!r}: before={b[f]!r} after={a[f]!r}"
            )
        assert b["_tags"] == a["_tags"], (
            f"record #{i} tag set diverged: before={b['_tags']!r} after={a['_tags']!r}"
        )


@pytest.mark.skipif(not os.path.exists(SRC), reason="smoke.bam fixture missing")
def test_field_level_header_references_preserved():
    """The header @SQ chain must survive the write — required for indexing."""
    bam_in = rubam.AlignmentFile(SRC, "rb")
    refs_before = list(bam_in.references)
    lens_before = list(bam_in.lengths)

    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "hdr.bam")
        bam_out = rubam.AlignmentFile(out, "wb", template=bam_in)
        for r in iter(bam_in):
            bam_out.write(r)
        bam_out.close()
        bam_in.close()

        bam_re = rubam.AlignmentFile(out, "rb")
        refs_after = list(bam_re.references)
        lens_after = list(bam_re.lengths)
        bam_re.close()

    assert refs_after == refs_before, (
        f"@SQ names diverged: before={refs_before} after={refs_after}"
    )
    assert lens_after == lens_before, (
        f"@SQ lengths diverged: before={lens_before} after={lens_after}"
    )
