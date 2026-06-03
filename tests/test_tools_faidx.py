import os
import tempfile
import rubam.tools

# Tiny FASTA we generate inline
FASTA_CONTENT = (
    ">chr1\n"
    "ACGTACGTACGTACGT\n"
    ">chr2\n"
    "TTTTTTTTAAAAAAAA\n"
)

def test_faidx_writes_fai_and_subsets():
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, "tiny.fa")
        with open(fa, "w") as f:
            f.write(FASTA_CONTENT)
        # Build the .fai
        result = rubam.tools.faidx(fa)
        assert result is None
        assert os.path.exists(fa + ".fai")
        # Pull a subsequence — chr1:1-4 is the first 4 bases ACGT (1-based inclusive)
        seq = rubam.tools.faidx(fa, region="chr1:1-4")
        assert isinstance(seq, str)
        assert seq.upper() == "ACGT"

def test_faidx_chr2_subseq():
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, "tiny.fa")
        with open(fa, "w") as f:
            f.write(FASTA_CONTENT)
        rubam.tools.faidx(fa)
        seq = rubam.tools.faidx(fa, region="chr2:5-12")
        assert seq.upper() == "TTTTAAAA"
