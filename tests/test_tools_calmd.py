import os
import tempfile
from pathlib import Path

import rubam
import rubam.tools

EXAMPLE_BAM = str(Path(__file__).parent / "example.bam")


def test_calmd_runs_and_emits_bam():
    with tempfile.TemporaryDirectory() as td:
        # Build a fasta that covers chr1 — content doesn't matter for v0.2
        # because we only count I + D in CIGAR (no mismatch reconstruction).
        fa = os.path.join(td, "ref.fa")
        with open(fa, "w") as f:
            f.write(">chr1\n")
            # 250 Mbp of N is unrealistic but we never read base-by-base.
            f.write("N" * 250_000_000 + "\n")
        rubam.tools.faidx(fa)
        out = os.path.join(td, "with_md.bam")
        rubam.tools.calmd(EXAMPLE_BAM, fa, output=out)
        assert os.path.exists(out) and os.path.getsize(out) > 0
        with rubam.AlignmentFile(out, "rb") as bam:
            n = sum(1 for _ in bam)
        assert n > 0


def test_calmd_inserts_NM_tag():
    """NM tag should be present on every record after calmd, even if some
    inputs already had it (replaced)."""
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, "ref.fa")
        with open(fa, "w") as f:
            f.write(">chr1\n")
            f.write("N" * 250_000_000 + "\n")
        rubam.tools.faidx(fa)
        out = os.path.join(td, "with_md.bam")
        rubam.tools.calmd(EXAMPLE_BAM, fa, output=out)
        with rubam.AlignmentFile(out, "rb") as bam:
            for r in bam:
                assert r.has_tag("NM"), f"record {r.query_name} missing NM"
                nm = r.get_tag("NM")
                assert isinstance(nm, int) and nm >= 0
                # For v0.2 NM = total I + total D across the CIGAR.
                ct = r.cigartuples or []
                # CIGAR ops: I=1, D=2.
                expected_nm = sum(L for op, L in ct if op in (1, 2))
                assert nm == expected_nm, (
                    f"NM mismatch on {r.query_name}: got {nm}, expected {expected_nm}"
                )
