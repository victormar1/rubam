import os, shutil, tempfile
from pathlib import Path
import rubam, rubam.tools

EXAMPLE_BAM = str(Path(__file__).parent / "example.bam")

def test_index_creates_bai():
    with tempfile.TemporaryDirectory() as td:
        bam = os.path.join(td, "x.bam")
        shutil.copy(EXAMPLE_BAM, bam)
        bai = bam + ".bai"
        if os.path.exists(bai):
            os.remove(bai)
        rubam.tools.index(bam)
        assert os.path.exists(bai) and os.path.getsize(bai) > 0
