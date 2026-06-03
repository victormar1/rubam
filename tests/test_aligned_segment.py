from pathlib import Path
import rubam

EXAMPLE_BAM = str(Path(__file__).parent / "example.bam")

def first_record():
    with rubam.AlignmentFile(EXAMPLE_BAM, "rb") as bam:
        return next(iter(bam))

def test_query_name_present():
    r = first_record()
    assert isinstance(r.query_name, str) and len(r.query_name) > 0

def test_reference_id_and_name():
    r = first_record()
    assert r.reference_id == 0
    assert r.reference_name == "chr1"

def test_reference_start_and_end():
    r = first_record()
    assert r.reference_start is not None and r.reference_start >= 0
    assert r.reference_end is not None and r.reference_end > r.reference_start

def test_template_length_is_int():
    r = first_record()
    assert isinstance(r.template_length, int)

def test_mapping_quality_in_range():
    r = first_record()
    assert 0 <= r.mapping_quality <= 255

def test_boolean_flags():
    r = first_record()
    for name in [
        "is_paired", "is_proper_pair", "is_unmapped", "is_mate_unmapped",
        "is_reverse", "is_mate_reverse", "is_read1", "is_read2",
        "is_secondary", "is_qcfail", "is_duplicate", "is_supplementary",
    ]:
        v = getattr(r, name)
        assert isinstance(v, bool), (name, type(v))

def test_example_bam_first_record_is_paired():
    r = first_record()
    assert r.is_paired is True
    assert r.is_secondary is False
    assert r.is_supplementary is False

def test_cigarstring():
    r = first_record()
    cs = r.cigarstring
    assert cs is None or (isinstance(cs, str) and any(c in cs for c in "MIDNSHP=X"))

def test_cigartuples():
    r = first_record()
    ct = r.cigartuples
    assert ct is None or all(
        isinstance(t, tuple) and len(t) == 2 and isinstance(t[0], int) and isinstance(t[1], int)
        for t in ct
    )

def test_cigar_consistency_when_present():
    r = first_record()
    cs, ct = r.cigarstring, r.cigartuples
    if cs is not None and ct is not None:
        # Number of CIGAR ops in the string equals number of tuples.
        n_ops_in_string = sum(1 for c in cs if c.isalpha() or c == "=")
        assert len(ct) == n_ops_in_string

def test_query_sequence_and_length():
    r = first_record()
    seq = r.query_sequence
    if seq is not None:
        assert isinstance(seq, str)
        assert all(c in "ACGTN" for c in seq.upper())
        assert r.query_length == len(seq)

def test_query_qualities_length_matches_sequence():
    r = first_record()
    quals = r.query_qualities
    seq = r.query_sequence
    if quals is not None and seq is not None:
        assert len(quals) == len(seq)
        assert all(isinstance(q, int) and 0 <= q <= 93 for q in quals)

def test_tags_returns_list_of_tuples():
    r = first_record()
    tags = r.tags
    assert isinstance(tags, list)
    for entry in tags:
        assert isinstance(entry, tuple) and len(entry) == 2
        name, value = entry
        assert isinstance(name, str) and len(name) == 2

def test_get_tag_returns_value_or_raises():
    r = first_record()
    if r.has_tag("NM"):
        nm = r.get_tag("NM")
        assert isinstance(nm, int)
    import pytest
    if not r.has_tag("XX"):
        with pytest.raises(KeyError):
            r.get_tag("XX")

def test_get_blocks_returns_list_of_intervals():
    r = first_record()
    blocks = r.get_blocks()
    assert isinstance(blocks, list)
    for s, e in blocks:
        assert isinstance(s, int) and isinstance(e, int) and s <= e

def test_get_reference_positions_count_matches_aligned_bases():
    r = first_record()
    pos = r.get_reference_positions()
    assert isinstance(pos, list)
    if r.cigartuples is not None:
        # Number of positions equals the bp count of M/=/X ops (op codes 0, 7, 8)
        aligned = sum(L for op, L in r.cigartuples if op in (0, 7, 8))
        assert len(pos) == aligned

def test_get_overlap_with_self_region():
    r = first_record()
    s, e = r.reference_start, r.reference_end
    if s is not None and e is not None:
        assert r.get_overlap(s, e) > 0
        assert r.get_overlap(e + 1000, e + 2000) == 0
