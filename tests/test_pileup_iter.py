from pathlib import Path
import rubam

EXAMPLE_BAM = str(Path(__file__).parent / "example.bam")
CHROM = "chr1"
START_0 = 999_999     # 0-based
END_0 = 1_000_020

def test_pileup_returns_columns_with_pos_and_depth():
    with rubam.AlignmentFile(EXAMPLE_BAM, "rb") as bam:
        cols = list(bam.pileup(CHROM, START_0, END_0))
        assert all(c.depth >= 0 for c in cols)
        assert all(c.reference_pos >= START_0 and c.reference_pos < END_0 for c in cols)
        assert all(c.reference_name == CHROM for c in cols)

def test_pileup_truncate_default_is_true():
    with rubam.AlignmentFile(EXAMPLE_BAM, "rb") as bam:
        cols = list(bam.pileup(CHROM, START_0, END_0))
        assert all(START_0 <= c.reference_pos < END_0 for c in cols)
