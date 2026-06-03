import os
import tempfile
from pathlib import Path

import rubam
import rubam.tools

EXAMPLE_BAM = str(Path(__file__).parent / "example.bam")

def test_view_count_only_returns_int():
    n = rubam.tools.view(EXAMPLE_BAM, count_only=True)
    assert isinstance(n, int) and n > 0

def test_view_with_region_filter():
    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, "subset.bam")
        rubam.tools.view(EXAMPLE_BAM, region="chr1:999990-1000010", output=out)
        with rubam.AlignmentFile(out, "rb") as bam:
            n = sum(1 for _ in bam)
        assert n > 0

def test_view_min_mapq_filter():
    n_all = rubam.tools.view(EXAMPLE_BAM, count_only=True)
    # min_mapq=255 keeps only records with MAPQ>=255 (typically none)
    n_high = rubam.tools.view(EXAMPLE_BAM, count_only=True, min_mapq=255)
    assert n_high <= n_all
