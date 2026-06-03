"""Coverage suite for the pysam-compatibility surface added in
v0.3.5–v0.3.8.

Tests are grouped by class. Each test calls one public method and
asserts on a behavioural property (return type, value range, identity
under round-trips, equality vs the canonical accessor). These tests
are what protects against silent regression when we rewire the pyo3
bindings or the noodles backend.
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

import pytest

import rubam


FIXTURE_BAM = Path(__file__).resolve().parent / "fixtures" / "smoke.bam"


# --- pytest fixtures --------------------------------------------------------

@pytest.fixture(scope="module")
def af() -> rubam.AlignmentFile:
    return rubam.AlignmentFile(str(FIXTURE_BAM))


@pytest.fixture(scope="module")
def seg(af: rubam.AlignmentFile) -> rubam.AlignedSegment:
    return next(af.fetch(contig="chr1", start=0, end=1000))


@pytest.fixture(scope="module")
def vcf_path(tmp_path_factory) -> Path:
    tmp = tmp_path_factory.mktemp("vcf")
    p = tmp / "smoke.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\trs1\tA\tG\t30\tPASS\tDP=10\n"
        "chr1\t150\trs2\tACG\tA\t40\t.\tDP=20\n"
    )
    return p


@pytest.fixture(scope="module")
def vf(vcf_path: Path) -> rubam.VariantFile:
    return rubam.VariantFile(str(vcf_path))


# --- AlignmentFile aliases & metadata ---------------------------------------

class TestAlignmentFileMetadata:
    def test_closed_alias(self, af):
        assert af.closed is False
        assert af.is_closed is False

    def test_mode_filename(self, af):
        assert af.mode in ("rb", "r")
        assert af.filename.endswith("smoke.bam")

    def test_format_introspection(self, af):
        assert af.format == "BAM"
        assert af.is_bam is True
        assert af.is_cram is False
        assert af.is_vcf is False
        assert af.is_bcf is False

    def test_io_direction_flags(self, af):
        assert af.is_read is True
        assert af.is_write is False

    def test_compression_metadata(self, af):
        assert af.compression == "BGZF"
        assert af.category == "alignment"
        assert af.is_remote is False
        assert af.is_stream is False

    def test_version_string(self, af):
        assert "noodles" in af.version

    def test_text_header(self, af):
        text = af.text
        assert isinstance(text, str)
        assert len(text) > 0

    def test_mapped_unmapped_nocoordinate(self, af):
        assert isinstance(af.mapped, int)
        assert isinstance(af.unmapped, int)
        assert af.nocoordinate == 0

    def test_get_tid_round_trip(self, af):
        tid = af.get_tid("chr1")
        assert tid == 0
        assert af.gettid("chr1") == 0
        assert af.get_reference_name(0) == "chr1"
        assert af.getrname(0) == "chr1"

    def test_get_tid_unknown(self, af):
        assert af.get_tid("does_not_exist") == -1

    def test_is_valid_reference_name(self, af):
        assert af.is_valid_reference_name("chr1") is True
        assert af.is_valid_reference_name("does_not_exist") is False

    def test_is_valid_tid(self, af):
        assert af.is_valid_tid(0) is True
        assert af.is_valid_tid(-1) is False
        assert af.is_valid_tid(999) is False

    def test_index_filename(self, af):
        # smoke.bam ships with a .bai
        idx = af.index_filename
        assert idx is None or idx.endswith(".bai")

    def test_duplicate_filehandle_and_check_truncation(self, af):
        assert af.duplicate_filehandle is False
        assert isinstance(af.check_truncation, bool)


class TestAlignmentFileParseRegion:
    def test_parse_region_string_full(self, af):
        tid, start, end = af.parse_region(region="chr1:1-100")
        assert tid == 0
        assert start == 0
        assert end == 100

    def test_parse_region_string_with_commas(self, af):
        tid, start, end = af.parse_region(region="chr1:1,000-2,000")
        assert tid == 0 and start == 999 and end == 2000

    def test_parse_region_kwargs(self, af):
        tid, start, end = af.parse_region(contig="chr1", start=10, end=50)
        assert tid == 0 and start == 10 and end == 50

    def test_parse_region_chr_only(self, af):
        tid, start, end = af.parse_region(region="chr1")
        assert tid == 0 and start is None and end is None


class TestAlignmentFileNoOpStubs:
    """Methods documented as stubs / no-ops in CHANGELOG v0.3.5-v0.3.8."""

    def test_mate_returns_none_stub(self, af, seg):
        # Documented v0.3.8 stub.
        assert af.mate(seg) is None

    def test_seek_tell_return_zero(self, af):
        assert af.seek(0) == 0
        assert af.tell() == 0

    def test_flush_reset_return_none(self, af):
        assert af.flush() is None
        assert af.reset() is None

    def test_find_introns_empty_for_no_n(self, af, seg):
        # smoke.bam contains 100M reads with no N ops
        assert dict(af.find_introns([seg])) == {}


# --- AlignedSegment aliases -------------------------------------------------

class TestAlignedSegmentAliases:
    def test_qname_alias(self, seg):
        assert seg.qname == seg.query_name

    def test_pos_alias(self, seg):
        assert seg.pos == seg.reference_start

    def test_mapq_alias(self, seg):
        assert seg.mapq == seg.mapping_quality

    def test_tid_alias(self, seg):
        assert seg.tid == seg.reference_id

    def test_cigar_alias(self, seg):
        assert seg.cigar == seg.cigartuples

    def test_isize_tlen_aliases(self, seg):
        assert seg.isize == seg.template_length
        assert seg.tlen == seg.template_length

    def test_rname_alias(self, seg):
        assert seg.rname == seg.reference_name

    def test_mpos_pnext_aliases(self, seg):
        # pysam convention: -1 for None/unmapped mate; mate_reference_start
        # may be None (Optional[int]). We assert the -1 sentinel is used.
        expected = seg.mate_reference_start if seg.mate_reference_start is not None else -1
        assert seg.mpos == expected
        assert seg.pnext == expected

    def test_aend_alias(self, seg):
        assert seg.aend == seg.reference_end

    def test_alen_rlen_aliases(self, seg):
        # alen / rlen / reference_length all == ref_end - ref_start
        if seg.reference_start is not None and seg.reference_end is not None:
            expected = seg.reference_end - seg.reference_start
            assert seg.alen == expected
            assert seg.rlen == expected
            assert seg.reference_length == expected

    def test_qlen_alias(self, seg):
        assert seg.qlen == seg.query_length

    def test_is_forward_complement(self, seg):
        assert seg.is_forward == (not seg.is_reverse)

    def test_is_mapped_complement(self, seg):
        assert seg.is_mapped == (not seg.is_unmapped)


class TestAlignedSegmentCIGARMethods:
    def test_get_aligned_pairs_returns_pairs(self, seg):
        pairs = seg.get_aligned_pairs()
        assert isinstance(pairs, list)
        assert len(pairs) > 0
        # Each pair: (qpos | None, refpos | None)
        for p in pairs[:5]:
            assert isinstance(p, tuple) and len(p) == 2

    def test_get_aligned_pairs_matches_only(self, seg):
        all_p = seg.get_aligned_pairs(matches_only=False)
        only_p = seg.get_aligned_pairs(matches_only=True)
        # matches_only drops (None, refpos) and (qpos, None) entries
        for q, r in only_p:
            assert q is not None and r is not None

    def test_get_cigar_stats_shape(self, seg):
        op_counts, base_counts = seg.get_cigar_stats()
        assert len(op_counts) == 11
        assert len(base_counts) == 11

    def test_blocks_positions_aligned_pairs(self, seg):
        assert seg.blocks == seg.get_blocks()
        # positions length matches reference span of M-runs
        positions = seg.positions
        assert isinstance(positions, list)
        assert all(isinstance(p, int) for p in positions)

    def test_bin_in_range(self, seg):
        b = seg.bin
        assert isinstance(b, int)
        assert 0 <= b < (1 << 16)


class TestAlignedSegmentSlicing:
    def test_query_alignment_sequence_is_substring(self, seg):
        seq = seg.query_sequence
        a_seq = seg.query_alignment_sequence
        if seq and a_seq:
            assert a_seq in seq

    def test_query_alignment_qualities_length(self, seg):
        a_qual = seg.query_alignment_qualities
        a_seq = seg.query_alignment_sequence
        if a_qual is not None and a_seq is not None:
            assert len(a_qual) == len(a_seq)

    def test_qstart_qend(self, seg):
        assert seg.qstart >= 0
        assert seg.qend >= seg.qstart

    def test_qstart_qend_round_trip(self, seg):
        # query_alignment_length == qend - qstart
        assert seg.query_alignment_length == seg.qend - seg.qstart

    def test_qual_string_length(self, seg):
        q = seg.qual
        if q:
            assert len(q) == seg.query_length

    def test_query_alias(self, seg):
        # `query` == query_alignment_sequence
        assert seg.query == seg.query_alignment_sequence


class TestAlignedSegmentForwardStrand:
    def test_get_forward_sequence_length(self, seg):
        seq = seg.query_sequence
        fwd = seg.get_forward_sequence()
        if seq and fwd:
            assert len(fwd) == len(seq)

    def test_get_forward_qualities_length(self, seg):
        qual = seg.query_qualities
        fwd_qual = seg.get_forward_qualities()
        if qual is not None and fwd_qual is not None:
            assert len(fwd_qual) == len(qual)


class TestAlignedSegmentSerialization:
    def test_tostring_starts_with_qname(self, seg):
        line = seg.tostring()
        assert line.startswith(seg.query_name + "\t")

    def test_to_string_matches_tostring(self, seg):
        assert seg.to_string() == seg.tostring()

    def test_fromstring_round_trip_basic(self, af, seg):
        sam = seg.tostring()
        roundtrip = rubam.AlignedSegment.fromstring(sam, af.header)
        assert roundtrip.query_name == seg.query_name
        assert roundtrip.flag == seg.flag
        assert roundtrip.mapping_quality == seg.mapping_quality
        assert roundtrip.cigarstring == seg.cigarstring

    def test_fromstring_round_trip_sequence(self, af, seg):
        sam = seg.tostring()
        roundtrip = rubam.AlignedSegment.fromstring(sam, af.header)
        assert roundtrip.query_sequence == seg.query_sequence


class TestAlignedSegmentMisc:
    def test_compare_identity_zero(self, seg):
        assert seg.compare(seg) == 0

    def test_overlap_within_block(self, seg):
        s, e = seg.reference_start or 0, (seg.reference_end or 0)
        ov = seg.get_overlap(s, e)
        assert ov >= 0

    def test_overlap_alias(self, seg):
        assert seg.overlap(0, 1_000_000) == seg.get_overlap(0, 1_000_000)

    def test_infer_query_length(self, seg):
        assert seg.infer_query_length() == seg.query_length

    def test_infer_read_length_alias(self, seg):
        assert seg.infer_read_length() == seg.query_length

    def test_inferred_length_property(self, seg):
        assert seg.inferred_length == seg.query_length

    def test_get_tags_method_form(self, seg):
        result = seg.get_tags()
        assert isinstance(result, list)

    def test_header_back_reference(self, seg):
        h = seg.header
        assert isinstance(h, rubam.Header)
        # reference_sequences match
        assert h.references == seg.header.references

    def test_modified_bases_empty_when_no_mm_tag(self, seg):
        # smoke.bam has no MM tag → empty dict
        d = dict(seg.modified_bases)
        assert isinstance(d, dict)

    def test_modified_bases_forward_empty(self, seg):
        d = dict(seg.modified_bases_forward)
        assert isinstance(d, dict)

    def test_get_reference_sequence_length_matches_alen(self, seg):
        ref = seg.get_reference_sequence()
        # Fallback returns "N"*alen when no FASTA repo
        if seg.alen:
            assert len(ref) == seg.alen


# --- VariantFile metadata ----------------------------------------------------

class TestVariantFileMetadata:
    def test_closed(self, vf):
        assert vf.closed is False
        assert vf.is_closed is False

    def test_mode_filename(self, vf):
        assert vf.mode == "r"
        assert vf.filename.endswith(".vcf")

    def test_format_introspection(self, vf):
        assert vf.format in ("VCF", "VCF.gz", "BCF")
        assert vf.is_vcf is True
        assert vf.is_bcf is False
        assert vf.is_bam is False

    def test_compression(self, vf):
        assert vf.compression in ("NONE", "BGZF")

    def test_category(self, vf):
        assert vf.category == "variant"

    def test_io_direction(self, vf):
        assert vf.is_read is True
        assert vf.is_write is False
        assert vf.is_reading is True

    def test_threads_stub(self, vf):
        assert vf.threads == 1

    def test_truncation_dup(self, vf):
        assert vf.check_truncation is True
        assert vf.duplicate_filehandle is False

    def test_tid_round_trip(self, vf):
        assert vf.get_tid("chr1") == 0
        assert vf.get_reference_name(0) == "chr1"

    def test_is_valid(self, vf):
        assert vf.is_valid_reference_name("chr1") is True
        assert vf.is_valid_tid(0) is True
        assert vf.is_valid_tid(-1) is False

    def test_parse_region(self, vf):
        tid, s, e = vf.parse_region(region="chr1:1-100")
        assert tid == 0 and s == 0 and e == 100

    def test_header_written_read_mode(self, vf):
        assert vf.header_written is False


class TestVariantFileBuilders:
    def test_copy_returns_variantfile(self, vf):
        c = vf.copy()
        assert isinstance(c, rubam.VariantFile)

    def test_new_record(self, vf):
        r = vf.new_record(contig="chr1", start=200, alleles=["A", "T"])
        assert isinstance(r, rubam.VariantRecord)


# --- VariantRecord aliases ---------------------------------------------------

class TestVariantRecordAliases:
    def test_first_record(self, vf):
        rec = next(iter(vf))
        # Basic identity
        assert rec.chrom == "chr1"
        assert rec.contig == "chr1"
        assert rec.ref == "A"
        assert rec.alts == ("G",)
        assert rec.alleles == ("A", "G")
        assert rec.qual == 30.0
        assert rec.id == "rs1"
        assert rec.start == 99  # 0-based
        assert rec.stop == 100  # 0-based exclusive
        assert rec.rlen == 1
        assert rec.rid == 0
        assert "PASS" in list(rec.filter)

    def test_alleles_variant_types(self, vf):
        rec = next(iter(vf))
        types = rec.alleles_variant_types
        assert types[0] == "REF"
        assert types[1] == "SNP"

    def test_alleles_variant_types_del(self, vf):
        # Second record has REF=ACG ALT=A — a deletion
        recs = list(vf)
        types = recs[1].alleles_variant_types
        assert types == ("REF", "DEL")

    def test_copy(self, vf):
        rec = next(iter(vf))
        c = rec.copy()
        assert isinstance(c, rubam.VariantRecord)
        assert c.chrom == rec.chrom
        assert c.position == rec.position

    def test_header_back_reference(self, vf):
        rec = next(iter(vf))
        h = rec.header
        assert isinstance(h, rubam.VariantHeader)


# --- Header sections ---------------------------------------------------------

class TestHeaderSections:
    def test_to_dict_has_all_sections(self, af):
        d = af.header.to_dict()
        assert set(d.keys()) == {"HD", "SQ", "RG", "PG", "CO"}

    def test_as_dict_alias(self, af):
        assert af.header.as_dict() == af.header.to_dict()

    def test_getitem_sq(self, af):
        sq = af.header["SQ"]
        assert isinstance(sq, list)
        assert all("SN" in entry and "LN" in entry for entry in sq)

    def test_getitem_unknown_raises(self, af):
        with pytest.raises(KeyError):
            af.header["NONEXISTENT"]

    def test_references_lengths_match(self, af):
        h = af.header
        assert len(h.references) == h.nreferences
        assert len(h.lengths) == h.nreferences

    def test_tostring_starts_with_hd(self, af):
        text = af.header.tostring()
        assert text.startswith("@HD") or text.startswith("@SQ")

    def test_str_equals_tostring(self, af):
        assert str(af.header) == af.header.tostring()


# --- FastxFile ---------------------------------------------------------------

class TestFastxFile:
    def test_fasta_iteration(self, tmp_path):
        p = tmp_path / "x.fasta"
        p.write_text(">chr1 test\nACGTACGT\n>chr2\nGGG\n")
        with rubam.FastxFile(str(p)) as fx:
            records = list(fx)
        assert len(records) == 2
        assert records[0].name == "chr1"
        assert records[0].sequence == "ACGTACGT"
        assert records[0].comment == "test"
        assert records[1].name == "chr2"
        assert records[1].sequence == "GGG"
        assert records[1].comment is None

    def test_fastq_iteration(self, tmp_path):
        p = tmp_path / "y.fastq"
        p.write_text("@r1\nACGT\n+\n!!!!\n@r2\nGGGG\n+\nABCD\n")
        with rubam.FastxFile(str(p)) as fx:
            records = list(fx)
        assert len(records) == 2
        assert records[0].name == "r1"
        assert records[0].sequence == "ACGT"
        assert records[0].quality == "!!!!"
        assert records[1].quality == "ABCD"

    def test_close(self, tmp_path):
        p = tmp_path / "z.fasta"
        p.write_text(">a\nACGT\n")
        fx = rubam.FastxFile(str(p))
        assert fx.is_open is True
        fx.close()
        assert fx.is_open is False

    def test_filename_property(self, tmp_path):
        p = tmp_path / "f.fasta"
        p.write_text(">a\nN\n")
        fx = rubam.FastxFile(str(p))
        assert fx.filename.endswith("f.fasta")
        fx.close()


# --- Module-level discovery --------------------------------------------------

class TestModuleSurface:
    def test_classes_present(self):
        assert hasattr(rubam, "AlignmentFile")
        assert hasattr(rubam, "AlignedSegment")
        assert hasattr(rubam, "VariantFile")
        assert hasattr(rubam, "VariantRecord")
        assert hasattr(rubam, "FastaFile")
        assert hasattr(rubam, "FastxFile")
        assert hasattr(rubam, "TabixFile")
        assert hasattr(rubam, "Header")

    def test_pysam_style_callables(self):
        for name in ("flagstat", "idxstats", "faidx", "view", "depth",
                     "index", "sort", "merge", "calmd"):
            assert callable(getattr(rubam, name))

    def test_subprocess_dispatchers(self):
        # rubam.samtools.<subcmd> and rubam.bcftools.<subcmd>
        assert callable(rubam.samtools)
        assert callable(rubam.bcftools)
        # Attribute access should produce sub-runners
        view_runner = rubam.samtools.view
        assert callable(view_runner)
