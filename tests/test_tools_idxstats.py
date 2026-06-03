from pathlib import Path
import rubam.tools

EXAMPLE_BAM = str(Path(__file__).parent / "example.bam")

def test_idxstats_returns_per_chrom_rows():
    rows = rubam.tools.idxstats(EXAMPLE_BAM)
    assert isinstance(rows, list) and len(rows) >= 1
    for row in rows:
        assert {"contig", "length", "mapped", "unmapped"}.issubset(row.keys())
    chr1 = next(r for r in rows if r["contig"] == "chr1")
    assert chr1["length"] > 0
    assert chr1["mapped"] >= 0
