import rubam
import rubam.tools


def test_tools_module_imports():
    assert rubam.tools is not None


def test_tools_flagstat_forwards_to_v01():
    """rubam.tools.flagstat should produce the same dict as rubam.flag_stats."""
    from pathlib import Path
    bam = str(Path(__file__).parent / "example.bam")
    a = rubam.tools.flagstat(bam)
    b = rubam.flag_stats(bam)
    assert a == b
    assert "total" in a and a["total"] > 0
