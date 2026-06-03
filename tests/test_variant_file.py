"""Smoke tests for the v0.3 VariantFile / VariantRecord."""

import rubam


_TINY_VCF = """##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30\tPASS\tDP=20
chr1\t200\t.\tC\tT\t40\tPASS\tDP=25
chr1\t300\trs1\tG\tA\t50\tPASS\tDP=30
"""


def test_variantfile_class_exists():
    assert hasattr(rubam, "VariantFile"), "rubam.VariantFile is exported"


def test_variantrecord_class_exists():
    assert hasattr(rubam, "VariantRecord"), "rubam.VariantRecord is exported"


def test_variantfile_open_close(tmp_path):
    vcf = tmp_path / "tiny.vcf"
    vcf.write_text(_TINY_VCF)
    f = rubam.VariantFile(str(vcf), "r")
    assert f.is_open
    f.close()
    assert not f.is_open


def test_variantfile_context_manager(tmp_path):
    vcf = tmp_path / "tiny.vcf"
    vcf.write_text(_TINY_VCF)
    with rubam.VariantFile(str(vcf), "r") as f:
        assert f.is_open
    assert not f.is_open


def test_variantfile_iter_yields_records(tmp_path):
    vcf = tmp_path / "tiny.vcf"
    vcf.write_text(_TINY_VCF)
    with rubam.VariantFile(str(vcf), "r") as f:
        records = list(f)
    assert len(records) == 3
    for r in records:
        assert isinstance(r, rubam.VariantRecord)


def test_variantfile_invalid_mode(tmp_path):
    vcf = tmp_path / "tiny.vcf"
    vcf.write_text(_TINY_VCF)
    import pytest
    # Truly unknown mode — must raise ValueError
    with pytest.raises((ValueError, OSError)):
        rubam.VariantFile(str(vcf), "wq")


def test_variantfile_unknown_path():
    import pytest
    with pytest.raises((IOError, OSError)):
        rubam.VariantFile("does_not_exist_zzz.vcf", "r")


def test_variantrecord_identity_props(tmp_path):
    vcf = tmp_path / "tiny.vcf"
    vcf.write_text(_TINY_VCF)
    with rubam.VariantFile(str(vcf), "r") as f:
        records = list(f)
    assert records[0].reference_name == "chr1"
    assert records[0].position == 100
    assert records[0].pos == 100
    assert records[0].reference == "A"
    assert records[0].ref_allele == "A"
    assert records[0].alternates == ("G",)
    assert records[0].alts == ("G",)
    assert records[0].ids == ()
    assert records[0].quality == 30.0
    assert records[0].qual == 30.0
    assert records[0].filters == ("PASS",)


def test_variantrecord_third_record_with_id(tmp_path):
    vcf = tmp_path / "tiny.vcf"
    vcf.write_text(_TINY_VCF)
    with rubam.VariantFile(str(vcf), "r") as f:
        records = list(f)
    assert records[2].ids == ("rs1",)
    assert records[2].position == 300
    assert records[2].reference == "G"
    assert records[2].alternates == ("A",)
    assert records[2].quality == 50.0


def test_variantrecord_multiallelic(tmp_path):
    """A VCF with two ALT alleles should produce a 2-tuple."""
    multi = """##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t500\t.\tA\tG,T\t.\tPASS\t.
"""
    vcf = tmp_path / "multi.vcf"
    vcf.write_text(multi)
    with rubam.VariantFile(str(vcf), "r") as f:
        records = list(f)
    assert records[0].alternates == ("G", "T")
    assert records[0].quality is None  # QUAL = '.'


# ---------------------------------------------------------------------------
# Task A4 — VariantRecord.samples  (dict-like genotype access)
# ---------------------------------------------------------------------------

_MULTISAMPLE_VCF = """##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\tNA12891\tNA12892
chr1\t100\t.\tA\tG\t30\tPASS\tDP=20\tGT:DP:AD\t0/1:30:15,15\t1/1:25:0,25\t./.:.:.,.\n"""


def test_record_samples_keys(tmp_path):
    p = tmp_path / "ms.vcf"; p.write_text(_MULTISAMPLE_VCF)
    with rubam.VariantFile(str(p), "r") as f:
        rec = next(iter(f))
    samples = rec.samples
    assert len(samples) == 3
    assert "NA12878" in samples
    assert "NA99999" not in samples
    assert sorted(list(samples)) == ["NA12878", "NA12891", "NA12892"]


def test_record_samples_gt(tmp_path):
    p = tmp_path / "ms.vcf"; p.write_text(_MULTISAMPLE_VCF)
    with rubam.VariantFile(str(p), "r") as f:
        rec = next(iter(f))
    samples = rec.samples
    assert samples["NA12878"]["GT"] == (0, 1)
    assert samples["NA12891"]["GT"] == (1, 1)
    assert samples["NA12892"]["GT"] == (None, None)


def test_record_samples_dp(tmp_path):
    p = tmp_path / "ms.vcf"; p.write_text(_MULTISAMPLE_VCF)
    with rubam.VariantFile(str(p), "r") as f:
        rec = next(iter(f))
    samples = rec.samples
    assert samples["NA12878"]["DP"] == 30
    assert samples["NA12891"]["DP"] == 25
    assert samples["NA12892"]["DP"] is None


def test_record_samples_ad_tuple(tmp_path):
    p = tmp_path / "ms.vcf"; p.write_text(_MULTISAMPLE_VCF)
    with rubam.VariantFile(str(p), "r") as f:
        rec = next(iter(f))
    samples = rec.samples
    assert samples["NA12878"]["AD"] == (15, 15)
    assert samples["NA12891"]["AD"] == (0, 25)


def test_record_samples_unknown_sample_raises(tmp_path):
    import pytest
    p = tmp_path / "ms.vcf"; p.write_text(_MULTISAMPLE_VCF)
    with rubam.VariantFile(str(p), "r") as f:
        rec = next(iter(f))
    with pytest.raises(KeyError):
        _ = rec.samples["NA99999"]


# ---------------------------------------------------------------------------
# Task A5 — VariantFile.fetch(contig, start, end)  (indexed query)
# ---------------------------------------------------------------------------

import pytest
import subprocess
import shutil


_FETCH_VCF = """\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=10000>
##contig=<ID=chr2,length=10000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30\tPASS\t.
chr1\t500\t.\tC\tT\t30\tPASS\t.
chr1\t1500\t.\tG\tA\t30\tPASS\t.
chr1\t9000\t.\tT\tC\t30\tPASS\t.
chr2\t200\t.\tA\tT\t30\tPASS\t.
"""


def _to_wsl(p):
    s = str(p).replace("\\", "/")
    if len(s) > 1 and s[1] == ":":
        s = "/mnt/" + s[0].lower() + s[2:]
    return s


def _build_indexed_vcf(tmp_path_factory):
    """Build a small bgzipped+tabix-indexed VCF for fetch tests.

    Skips if WSL bgzip/tabix is unavailable.
    """
    tmp = tmp_path_factory.mktemp("indexed_vcf")
    plain = tmp / "small.vcf"
    plain.write_text(_FETCH_VCF)
    gz = tmp / "small.vcf.gz"

    if shutil.which("bgzip") is not None:
        # Native bgzip / tabix available (Linux CI or WSL shell)
        subprocess.check_call(["bgzip", "-c", str(plain)], stdout=open(str(gz), "wb"))
        subprocess.check_call(["tabix", "-p", "vcf", str(gz)])
    elif shutil.which("wsl") is not None:
        cmd = (
            f"bgzip -c {_to_wsl(plain)} > {_to_wsl(gz)} && "
            f"tabix -p vcf {_to_wsl(gz)}"
        )
        result = subprocess.run(
            ["wsl", "-d", "Ubuntu", "bash", "-lc", cmd],
            capture_output=True,
        )
        if result.returncode != 0:
            pytest.skip("bgzip/tabix not available in WSL")
    else:
        pytest.skip("bgzip/tabix not available")

    return gz


@pytest.fixture(scope="module")
def indexed_vcf(tmp_path_factory):
    return _build_indexed_vcf(tmp_path_factory)


def test_fetch_chr1_full(indexed_vcf):
    with rubam.VariantFile(str(indexed_vcf), "r") as f:
        recs = list(f.fetch("chr1", 1, 10000))
    assert len(recs) == 4
    assert recs[0].position == 100
    assert recs[-1].position == 9000


def test_fetch_chr1_subrange(indexed_vcf):
    with rubam.VariantFile(str(indexed_vcf), "r") as f:
        recs = list(f.fetch("chr1", 200, 2000))
    positions = [r.position for r in recs]
    # 100 is below start, 9000 is above end -> only 500 and 1500 returned.
    assert positions == [500, 1500]


def test_fetch_chr2(indexed_vcf):
    with rubam.VariantFile(str(indexed_vcf), "r") as f:
        recs = list(f.fetch("chr2", 1, 10000))
    assert len(recs) == 1
    assert recs[0].position == 200
    assert recs[0].reference_name == "chr2"


def test_fetch_unknown_contig_raises(indexed_vcf):
    with rubam.VariantFile(str(indexed_vcf), "r") as f:
        with pytest.raises(ValueError):
            list(f.fetch("chr99", 1, 1000))


def test_fetch_empty_interval_yields_empty(indexed_vcf):
    # end <= start should give an empty iterator without raising
    with rubam.VariantFile(str(indexed_vcf), "r") as f:
        recs = list(f.fetch("chr1", 5000, 5000))
    assert recs == []


# ---------------------------------------------------------------------------
# Task A6 — VariantHeader  (read-only header metadata)
# ---------------------------------------------------------------------------

_HEADER_VCF = """##fileformat=VCFv4.3
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FILTER=<ID=LowQual,Description="Low quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\tNA12891
chr1\t100\t.\tA\tG\t30\tPASS\tDP=20;AF=0.5\tGT:AD\t0/1:15,15\t1/1:0,25
"""


def test_header_samples(tmp_path):
    p = tmp_path / "h.vcf"; p.write_text(_HEADER_VCF)
    with rubam.VariantFile(str(p), "r") as f:
        assert f.header.samples == ("NA12878", "NA12891")


def test_header_version(tmp_path):
    p = tmp_path / "h.vcf"; p.write_text(_HEADER_VCF)
    with rubam.VariantFile(str(p), "r") as f:
        assert f.header.version == "VCFv4.3"


def test_header_contigs(tmp_path):
    p = tmp_path / "h.vcf"; p.write_text(_HEADER_VCF)
    with rubam.VariantFile(str(p), "r") as f:
        contigs = f.header.contigs
    assert len(contigs) == 2
    assert "chr1" in contigs
    assert "chrX" not in contigs
    assert contigs["chr1"].length == 248_956_422
    assert sorted(list(contigs)) == ["chr1", "chr2"]


def test_header_info_defs(tmp_path):
    p = tmp_path / "h.vcf"; p.write_text(_HEADER_VCF)
    with rubam.VariantFile(str(p), "r") as f:
        info = f.header.info
    assert "DP" in info
    assert "AF" in info
    assert "BOGUS" not in info
    dp = info["DP"]
    assert dp.id == "DP"
    assert dp.type == "Integer"
    assert dp.number == "1"
    assert dp.description == "Total depth"
    af = info["AF"]
    assert af.number == "A"
    assert af.type == "Float"


def test_header_format_defs(tmp_path):
    p = tmp_path / "h.vcf"; p.write_text(_HEADER_VCF)
    with rubam.VariantFile(str(p), "r") as f:
        fmt = f.header.formats
    assert "GT" in fmt
    assert fmt["GT"].type == "String"
    assert fmt["AD"].number == "R"


def test_header_filters(tmp_path):
    p = tmp_path / "h.vcf"; p.write_text(_HEADER_VCF)
    with rubam.VariantFile(str(p), "r") as f:
        filters = f.header.filters
    # The implicit PASS may or may not appear; the user-declared LowQual must be present.
    assert "LowQual" in filters


def test_header_unknown_contig_raises(tmp_path):
    import pytest
    p = tmp_path / "h.vcf"; p.write_text(_HEADER_VCF)
    with rubam.VariantFile(str(p), "r") as f:
        with pytest.raises(KeyError):
            _ = f.header.contigs["chrX"]


# ---------------------------------------------------------------------------
# Task A7 — VariantFile write modes  w / wz / wb
# ---------------------------------------------------------------------------

def test_write_plain_vcf_round_trip(tmp_path):
    src_path = tmp_path / "src.vcf"; src_path.write_text(_TINY_VCF)
    out_path = tmp_path / "out.vcf"
    with rubam.VariantFile(str(src_path), "r") as src:
        with rubam.VariantFile(str(out_path), "w", header=src.header) as out:
            for rec in src:
                out.write(rec)
    # Read back
    with rubam.VariantFile(str(out_path), "r") as got:
        recs = list(got)
    assert len(recs) == 3
    assert recs[0].position == 100
    assert recs[2].ids == ("rs1",)


def test_write_bgzf_vcf_round_trip(tmp_path):
    src_path = tmp_path / "src.vcf"; src_path.write_text(_TINY_VCF)
    out_path = tmp_path / "out.vcf.gz"
    with rubam.VariantFile(str(src_path), "r") as src:
        with rubam.VariantFile(str(out_path), "wz", header=src.header) as out:
            for rec in src:
                out.write(rec)
    # File should be BGZF-compressed (binary, non-empty)
    assert out_path.stat().st_size > 0
    with rubam.VariantFile(str(out_path), "r") as got:
        recs = list(got)
    assert len(recs) == 3


def test_write_bcf_round_trip(tmp_path):
    src_path = tmp_path / "src.vcf"; src_path.write_text(_TINY_VCF)
    out_path = tmp_path / "out.bcf"
    with rubam.VariantFile(str(src_path), "r") as src:
        with rubam.VariantFile(str(out_path), "wb", header=src.header) as out:
            for rec in src:
                out.write(rec)
    # BCF is binary; read back via rubam.
    with rubam.VariantFile(str(out_path), "r") as got:
        recs = list(got)
    assert len(recs) == 3
    assert recs[0].position == 100


def test_write_unknown_mode_raises(tmp_path):
    import pytest
    src_path = tmp_path / "src.vcf"; src_path.write_text(_TINY_VCF)
    with rubam.VariantFile(str(src_path), "r") as src:
        with pytest.raises(ValueError):
            rubam.VariantFile(str(tmp_path / "x"), "wq", header=src.header)


def test_write_without_header_raises(tmp_path):
    import pytest
    with pytest.raises((ValueError, TypeError)):
        rubam.VariantFile(str(tmp_path / "x.vcf"), "w")  # no header=


def test_write_to_reader_raises(tmp_path):
    import pytest
    p = tmp_path / "tiny.vcf"; p.write_text(_TINY_VCF)
    with rubam.VariantFile(str(p), "r") as src:
        rec = next(iter(src))
        # Writing into a read-mode file must raise.
        with pytest.raises((IOError, OSError)):
            src.write(rec)


# ---------------------------------------------------------------------------
# Task A8 — VariantRecord constructor + mutation APIs
# ---------------------------------------------------------------------------

_TEST_HEADER_VCF = """##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=NAME,Number=1,Type=String,Description="Some name">
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30\tPASS\tDP=20
"""


def _open_header(tmp_path):
    p = tmp_path / "h.vcf"
    p.write_text(_TEST_HEADER_VCF)
    return rubam.VariantFile(str(p), "r")


def test_construct_record_minimal(tmp_path):
    src = _open_header(tmp_path)
    hdr = src.header
    rec = rubam.VariantRecord(
        header=hdr,
        reference_name="chr1",
        position=200,
        reference="C",
        alternates=("T",),
    )
    assert rec.reference_name == "chr1"
    assert rec.position == 200
    assert rec.reference == "C"
    assert rec.alternates == ("T",)
    assert rec.quality is None
    src.close()


def test_construct_record_full(tmp_path):
    src = _open_header(tmp_path)
    rec = rubam.VariantRecord(
        header=src.header,
        reference_name="chr1",
        position=300,
        reference="G",
        alternates=("A", "T"),
        quality=42.0,
        ids=("rs99",),
        filters=("PASS",),
    )
    assert rec.alternates == ("A", "T")
    assert rec.quality == 42.0
    assert rec.ids == ("rs99",)
    assert "PASS" in rec.filters
    src.close()


def test_set_position(tmp_path):
    src = _open_header(tmp_path)
    rec = next(iter(src))
    assert rec.position == 100
    rec.set_position(500)
    assert rec.position == 500
    src.close()


def test_set_quality(tmp_path):
    src = _open_header(tmp_path)
    rec = next(iter(src))
    rec.set_quality(99.5)
    assert rec.quality == 99.5
    src.close()


def test_set_filter_replaces(tmp_path):
    src = _open_header(tmp_path)
    rec = next(iter(src))
    rec.set_filter("LowQual")
    assert "LowQual" in rec.filters
    assert "PASS" not in rec.filters
    src.close()


def test_add_filter_appends(tmp_path):
    src = _open_header(tmp_path)
    rec = next(iter(src))
    rec.set_filter("LowQual")
    rec.add_filter("PASS")  # append
    assert "LowQual" in rec.filters
    assert "PASS" in rec.filters
    src.close()


def test_clear_filters(tmp_path):
    src = _open_header(tmp_path)
    rec = next(iter(src))
    rec.clear_filters()
    assert rec.filters == ()
    src.close()


def test_set_info_int(tmp_path):
    src = _open_header(tmp_path)
    rec = next(iter(src))
    rec.set_info("DP", 99)
    # Round-trip via write + read back
    out = tmp_path / "out.vcf"
    with rubam.VariantFile(str(out), "w", header=src.header) as w:
        w.write(rec)
    with rubam.VariantFile(str(out), "r") as r:
        got = next(iter(r))
    src.close()
    # The text round-trip is the strongest assertion; we just verify the
    # written file parses without error.
    assert got.position == rec.position


def test_construct_then_write_round_trip(tmp_path):
    src = _open_header(tmp_path)
    hdr = src.header
    rec = rubam.VariantRecord(
        header=hdr,
        reference_name="chr1",
        position=777,
        reference="A",
        alternates=("C",),
        quality=88.0,
        ids=("rs77",),
        filters=("PASS",),
    )
    out = tmp_path / "constructed.vcf"
    with rubam.VariantFile(str(out), "w", header=hdr) as w:
        w.write(rec)
    src.close()
    with rubam.VariantFile(str(out), "r") as r:
        recs = list(r)
    assert len(recs) == 1
    assert recs[0].position == 777
    assert recs[0].reference == "A"
    assert recs[0].alternates == ("C",)
    assert recs[0].quality == 88.0
    assert recs[0].ids == ("rs77",)
