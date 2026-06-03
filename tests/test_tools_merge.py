import os, shutil, tempfile
from pathlib import Path

import rubam
import rubam.tools

EXAMPLE_BAM = str(Path(__file__).parent / "example.bam")

def test_merge_two_copies_doubles_record_count():
    with tempfile.TemporaryDirectory() as td:
        a, b = os.path.join(td, "a.bam"), os.path.join(td, "b.bam")
        shutil.copy(EXAMPLE_BAM, a)
        shutil.copy(EXAMPLE_BAM, b)
        out = os.path.join(td, "merged.bam")
        rubam.tools.merge([a, b], out)
        with rubam.AlignmentFile(EXAMPLE_BAM, "rb") as src:
            n_src = sum(1 for _ in src)
        with rubam.AlignmentFile(out, "rb") as dst:
            n_dst = sum(1 for _ in dst)
        assert n_dst == 2 * n_src

def test_merge_zero_inputs_raises():
    import pytest
    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, "out.bam")
        with pytest.raises((IOError, OSError, ValueError)):
            rubam.tools.merge([], out)
