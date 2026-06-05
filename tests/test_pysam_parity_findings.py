"""Regression tests for three pysam-parity findings on real-world hg38 BAMs.

1. Header parser must accept real-world hg38 htslib/GATK headers that
   noodles' strict parser rejects: @HD without VN, multi-part versions
   (VN:1.6.0), duplicate @PG/@RG IDs.
2. count_coverage must match pysam: quality_threshold default 15 (>= test),
   read_callback='all' (0x704 mask), and no depth cap.
3. count must match pysam: read_callback default 'nofilter' (counts secondary,
   supplementary, duplicate, QC-fail), 'all' applies the 0x704 mask.

Expected count/coverage values were captured from pysam 0.24.0 on the same
fixture (tests/fixtures/pysam_parity.bam).
"""
import struct
import os
import tempfile

import rubam

HERE = os.path.dirname(os.path.abspath(__file__))
PARITY_BAM = os.path.join(HERE, "fixtures", "pysam_parity.bam")


# --------------------------------------------------------------------------
# Finding 1 — tolerant header parsing
# --------------------------------------------------------------------------
def _build_raw_bam(sam_text, refs):
    out = bytearray(b"BAM\x01")
    text = sam_text.encode()
    out += struct.pack("<i", len(text)) + text
    out += struct.pack("<i", len(refs))
    for name, ln in refs:
        nm = name.encode() + b"\x00"
        out += struct.pack("<i", len(nm)) + nm + struct.pack("<i", ln)
    return bytes(out)


def _write_header_only_bam(sam_text, refs):
    path = os.path.join(tempfile.gettempdir(), "rubam_hdr_parity.bam")
    raw = _build_raw_bam(sam_text, refs)
    f = rubam.BGZFile(path, "wb")
    f.write(raw)
    f.close()
    return path


# A representative GRCh38-style contig subset (names/lengths real).
REFS = [
    ("chr1", 248956422),
    ("chr2", 242193529),
    ("chrM", 16569),
    ("chrUn_KI270742v1", 186739),
]
SQ = "".join(f"@SQ\tSN:{n}\tLN:{l}\n" for n, l in REFS)


def _assert_opens(sam_text):
    path = _write_header_only_bam(sam_text, REFS)
    with rubam.AlignmentFile(path, "rb") as af:
        assert af.nreferences == len(REFS)
        assert tuple(af.references) == tuple(n for n, _ in REFS)
        assert tuple(af.lengths) == tuple(l for _, l in REFS)


def test_header_hd_without_vn():
    _assert_opens("@HD\tSO:coordinate\n" + SQ)


def test_header_multipart_version():
    _assert_opens("@HD\tVN:1.6.0\tSO:coordinate\n" + SQ)


def test_header_no_hd_line():
    _assert_opens(SQ)


def test_header_duplicate_pg_id():
    txt = ("@HD\tVN:1.6\n" + SQ +
           "@PG\tID:bwa\tPN:bwa\tCL:bwa mem ref.fa r1.fq r2.fq\n"
           "@PG\tID:bwa\tPN:bwa\tCL:bwa mem ref.fa (rerun)\n")
    _assert_opens(txt)


def test_header_well_formed_still_works():
    txt = ("@HD\tVN:1.6\tSO:coordinate\n" + SQ +
           "@RG\tID:s1\tSM:NA12878\tPL:ILLUMINA\n"
           "@PG\tID:bwa\tPN:bwa\tVN:0.7.17\tCL:bwa mem\n")
    _assert_opens(txt)


# --------------------------------------------------------------------------
# Findings 2 & 3 — count / count_coverage parity with pysam
# (reference values from pysam 0.24.0 on tests/fixtures/pysam_parity.bam)
# --------------------------------------------------------------------------
def test_count_default_is_nofilter():
    # pysam count default read_callback='nofilter': every read in the region,
    # including secondary/supplementary/duplicate/qcfail.
    with rubam.AlignmentFile(PARITY_BAM, "rb") as af:
        assert af.count("chr1", 100, 120) == 10
        assert af.count("chr1") == 10


def test_count_read_callback_all():
    # 'all' skips UNMAP|SECONDARY|QCFAIL|DUP (0x704), keeps supplementary.
    with rubam.AlignmentFile(PARITY_BAM, "rb") as af:
        assert af.count("chr1", 100, 120, read_callback="all") == 7
        assert af.count("chr1", 100, 120, read_callback="nofilter") == 10


def test_free_count_reads_default_diverges_from_method():
    # The "count trap": the free count_reads() defaults to the samtools 0x704
    # mask (excludes secondary/dup/qcfail/unmap), whereas AlignmentFile.count()
    # defaults to pysam's 'nofilter' (counts everything). They diverge on the
    # same region by design. count_reads is 1-based inclusive; the method is
    # 0-based half-open, so count_reads(101, 120) covers af.count(.., 100, 120).
    free_default = rubam.count_reads(PARITY_BAM, "chr1", 101, 120)
    free_nofilter = rubam.count_reads(PARITY_BAM, "chr1", 101, 120, flag_filtered=0)
    with rubam.AlignmentFile(PARITY_BAM, "rb") as af:
        method_default = af.count("chr1", 100, 120)
    assert free_default == 7              # samtools 0x704 default
    assert method_default == 10           # pysam nofilter default
    assert free_default < method_default  # the documented trap
    # documented escape hatch: flag_filtered=0 reconciles the free fn with the method
    assert free_nofilter == method_default == 10


def test_count_coverage_default_quality_threshold_15():
    # First position has 5 reads with base-quality >= 15 (q15,q16,q20,q20 + the
    # boundary q15 read) under the default 'all' filter.
    with rubam.AlignmentFile(PARITY_BAM, "rb") as af:
        a, c, g, t = af.count_coverage("chr1", 100, 120)
        depth0 = a[0] + c[0] + g[0] + t[0]
        assert depth0 == 5


def test_count_coverage_quality_threshold_boundary():
    # >= semantics: qt=15 counts the q15 read (5), qt=16 drops it (4).
    with rubam.AlignmentFile(PARITY_BAM, "rb") as af:
        d15 = sum(arr[0] for arr in af.count_coverage("chr1", 100, 120, quality_threshold=15))
        d16 = sum(arr[0] for arr in af.count_coverage("chr1", 100, 120, quality_threshold=16))
        d0 = sum(arr[0] for arr in af.count_coverage("chr1", 100, 120, quality_threshold=0))
        assert d15 == 5
        assert d16 == 4
        assert d0 == 6  # 'all' filter: 6 mapped non-sec/dup/qcfail reads at pos 0


def test_count_coverage_read_callback_nofilter():
    with rubam.AlignmentFile(PARITY_BAM, "rb") as af:
        d = sum(arr[0] for arr in
                af.count_coverage("chr1", 100, 120, quality_threshold=0, read_callback="nofilter"))
        # nofilter at pos 0 adds secondary+duplicate+qcfail (3 more) -> 9
        assert d == 9
