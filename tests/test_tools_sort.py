import os
import shutil
import tempfile
from pathlib import Path

import rubam
import rubam.tools

EXAMPLE_BAM = str(Path(__file__).parent / "example.bam")

def test_sort_produces_coordinate_sorted_bam():
    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, "sorted.bam")
        rubam.tools.sort(EXAMPLE_BAM, out, threads=2)
        assert os.path.exists(out) and os.path.getsize(out) > 0
        with rubam.AlignmentFile(out, "rb") as bam:
            prev_rid, prev_pos = -1, -1
            for r in bam:
                rid = r.reference_id if r.reference_id is not None else -1
                pos = r.reference_start if r.reference_start is not None else -1
                if rid == prev_rid:
                    assert pos >= prev_pos, (rid, pos, prev_pos)
                else:
                    assert rid >= prev_rid
                prev_rid, prev_pos = rid, pos
