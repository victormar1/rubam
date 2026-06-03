from pathlib import Path
import os
import pytest
import rubam

EXAMPLE_BAM = str(Path(__file__).parent / "example.bam")

# Optional large CRAM + GRCh38 reference for local CRAM checks. Not bundled;
# point these env vars at local copies to exercise the CRAM path. Tests that
# need them skip when unset/missing.
_CRAM = os.environ.get("RUBAM_TEST_CRAM", "")
_REF = os.environ.get("RUBAM_TEST_REF", "")

def test_alignmentfile_class_exists():
    assert hasattr(rubam, "AlignmentFile"), "rubam.AlignmentFile is exported"

def test_alignedsegment_class_exists():
    assert hasattr(rubam, "AlignedSegment"), "rubam.AlignedSegment is exported"

def test_alignmentfile_open_close():
    bam = rubam.AlignmentFile(EXAMPLE_BAM, "rb")
    assert bam.is_open
    bam.close()
    assert not bam.is_open

def test_alignmentfile_context_manager():
    with rubam.AlignmentFile(EXAMPLE_BAM, "rb") as bam:
        assert bam.is_open
    assert not bam.is_open

def test_alignmentfile_references():
    with rubam.AlignmentFile(EXAMPLE_BAM, "rb") as bam:
        assert isinstance(bam.references, tuple)
        assert "chr1" in bam.references

def test_alignmentfile_lengths_matches_references():
    with rubam.AlignmentFile(EXAMPLE_BAM, "rb") as bam:
        assert isinstance(bam.lengths, tuple)
        assert len(bam.lengths) == len(bam.references)
        assert all(L > 0 for L in bam.lengths)

def test_alignmentfile_nreferences():
    with rubam.AlignmentFile(EXAMPLE_BAM, "rb") as bam:
        assert bam.nreferences == len(bam.references)

def test_alignmentfile_header_dict_round_trip():
    with rubam.AlignmentFile(EXAMPLE_BAM, "rb") as bam:
        h = bam.header
        d = h.to_dict()
        assert any(rec.get("SN") == "chr1" for rec in d.get("SQ", []))

def test_fetch_returns_iterator_of_records_in_region():
    with rubam.AlignmentFile(EXAMPLE_BAM, "rb") as bam:
        # example.bam has reads near chr1:1_000_000
        reads = list(bam.fetch("chr1", 999_990, 1_000_010))
        assert len(reads) > 0
        for r in reads:
            assert r.reference_name == "chr1"
            assert r.reference_start is not None
            # 0-based half-open overlap
            assert r.reference_end > 999_990
            assert r.reference_start < 1_000_010

def test_fetch_unknown_chromosome_raises():
    with rubam.AlignmentFile(EXAMPLE_BAM, "rb") as bam:
        with pytest.raises((ValueError, OSError)):
            list(bam.fetch("chrZZZ", 0, 100))

def test_fetch_empty_region_returns_no_records():
    with rubam.AlignmentFile(EXAMPLE_BAM, "rb") as bam:
        reads = list(bam.fetch("chr1", 1_000_000_000, 1_000_000_100))
        assert reads == []

def test_has_index():
    with rubam.AlignmentFile(EXAMPLE_BAM, "rb") as bam:
        assert bam.has_index() is True

def test_check_index_does_not_raise():
    with rubam.AlignmentFile(EXAMPLE_BAM, "rb") as bam:
        bam.check_index()  # returns None or True; no exception

def test_get_index_statistics_returns_per_chrom_rows():
    with rubam.AlignmentFile(EXAMPLE_BAM, "rb") as bam:
        rows = bam.get_index_statistics()
        assert isinstance(rows, list) and len(rows) >= 1
        for row in rows:
            assert "contig" in row and "mapped" in row and "unmapped" in row

def test_head_returns_first_n():
    with rubam.AlignmentFile(EXAMPLE_BAM, "rb") as bam:
        first5 = bam.head(5)
        assert isinstance(first5, list) and len(first5) <= 5
        for r in first5:
            assert isinstance(r, rubam.AlignedSegment)

def test_count_with_region():
    with rubam.AlignmentFile(EXAMPLE_BAM, "rb") as bam:
        n = bam.count("chr1", 999_990, 1_000_020)
        assert isinstance(n, int) and n > 0

def test_count_coverage_returns_4_lists():
    with rubam.AlignmentFile(EXAMPLE_BAM, "rb") as bam:
        a, c, g, t = bam.count_coverage("chr1", 1_000_000, 1_000_010)
        for arr in (a, c, g, t):
            assert isinstance(arr, list) and len(arr) == 10

def test_alignmentfile_is_open_works_for_unindexed_bam(tmp_path):
    """Reopening an unindexed BAM (e.g. fresh from rubam.tools.sort) must
    report is_open=True. Regression for the B2 streaming-fallback bug."""
    import shutil
    import rubam.tools

    src = str(Path(__file__).parent / "example.bam")
    out = str(tmp_path / "sorted_no_bai.bam")
    rubam.tools.sort(src, out)
    # No .bai is built; sort only sorts.
    assert not (Path(out + ".bai")).exists()
    bam = rubam.AlignmentFile(out, "rb")
    try:
        assert bam.is_open is True, "is_open should be True for an open unindexed BAM"
    finally:
        bam.close()
    assert bam.is_open is False, "is_open should be False after close()"


# ---------------------------------------------------------------------------
# CRAM smoke tests — skipped automatically if the local dataset is absent.
# ---------------------------------------------------------------------------

_CRAM_SKIP = pytest.mark.skipif(
    not os.path.exists(_CRAM) or not os.path.exists(_CRAM + ".crai"),
    reason="NA12878 CRAM dataset not available locally",
)
_REF_SKIP = pytest.mark.skipif(
    not os.path.exists(_REF) or not os.path.exists(_REF + ".fai"),
    reason="GRCh38 reference FASTA not indexed locally (run: samtools faidx)",
)


@_CRAM_SKIP
def test_cram_open_and_header():
    """CRAM: AlignmentFile opens successfully and returns a non-empty header."""
    with rubam.AlignmentFile(_CRAM) as f:
        assert f.is_open
        assert f.nreferences > 0
        assert isinstance(f.references, tuple) and len(f.references) > 0
        h = f.header
        d = h.to_dict()
        assert len(d.get("SQ", [])) > 0


@_CRAM_SKIP
def test_cram_has_index():
    """CRAM: has_index() returns True when .crai is present."""
    with rubam.AlignmentFile(_CRAM) as f:
        assert f.has_index() is True


@_CRAM_SKIP
@pytest.mark.xfail(
    reason=(
        "NA12878 NYGC 30x CRAM uses Huffman byte-series encoding not yet "
        "implemented in noodles-cram 0.90 (decode_take todo!()). "
        "AlignmentFile.open() and .header work; .fetch() panics inside noodles. "
        "Will become a real pass once noodles ships the codec."
    ),
    strict=True,
    raises=BaseException,
)
def test_cram_fetch_returns_records():
    """CRAM: fetch() API wiring + AlignedSegment type checks."""
    ref = _REF if os.path.exists(_REF) else None
    with rubam.AlignmentFile(_CRAM, reference_filename=ref) as f:
        refs = f.references
        assert f.nreferences > 0
        chrom = "chr1" if "chr1" in refs else "1"
        n = 0
        for r in f.fetch(chrom, 1_000_000, 1_100_000):
            assert isinstance(r, rubam.AlignedSegment)
            assert r.reference_name == chrom
            assert r.reference_start is not None
            n += 1
            if n >= 100:
                break
        assert n > 0, f"No reads returned for {chrom}:1_000_000-1_100_000"


@_CRAM_SKIP
@pytest.mark.xfail(
    reason=(
        "NA12878 NYGC 30x CRAM uses Huffman byte-series encoding not yet "
        "implemented in noodles-cram 0.90."
    ),
    strict=True,
    raises=BaseException,
)
def test_cram_aligned_segment_attributes():
    """CRAM: AlignedSegment attribute checks (xfail until noodles Huffman codec)."""
    ref = _REF if os.path.exists(_REF) else None
    with rubam.AlignmentFile(_CRAM, reference_filename=ref) as f:
        refs = f.references
        chrom = "chr1" if "chr1" in refs else "1"
        for r in f.fetch(chrom, 1_000_000, 1_100_000):
            assert isinstance(r.flag, int)
            assert isinstance(r.mapping_quality, int)
            assert r.query_name is not None
            assert r.reference_start is not None
            break
