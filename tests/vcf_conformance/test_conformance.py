"""VCF conformance suite: rubam.VariantFile vs pysam.VariantFile.

For each curated fixture under ``tests/vcf_conformance/fixtures/`` this module
compares the records returned by rubam against the records returned by
pysam, field by field. Known divergences (features that rubam does not yet
expose on the read-side) are marked with ``pytest.xfail`` -- never with a
silent ``skip`` -- so that the manuscript matrix stays honest.

Tested dimensions
-----------------
* core 7-tuple: (chrom, pos, ref, alts, ids, qual, filters)
* multi-allelic ALTs and Number=A/R/G INFO/FORMAT fields
* phased vs unphased GT
* missing values ("./.", ".")
* multi-sample FORMAT extraction
* complex FORMAT layouts (GT:DP:AD:GQ:PL)

NOT tested (see ``docs/vcf_conformance_matrix.md``): symbolic ALT (``<DEL>``),
breakends (BND), polyploid genotypes. Those are explicitly out of scope for
the current rubam v0.3 surface.

The suite is intentionally skipped wholesale when ``pysam`` is not
importable -- e.g. on Windows where pysam ships no wheel. The CI matrix runs
this file under WSL/Linux where pysam is available.
"""

from __future__ import annotations

import os
import pathlib

import pytest

pysam = pytest.importorskip(
    "pysam",
    reason="pysam not installed (expected on Windows; conformance runs on Linux/WSL)",
)

import rubam  # noqa: E402  -- after pysam guard


FIXTURE_DIR = pathlib.Path(__file__).parent / "fixtures"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _pysam_ids(rec) -> tuple:
    """Return the IDs field as a tuple, matching the rubam convention."""
    if rec.id is None:
        return ()
    # pysam returns a single string; multiple IDs are ";"-joined per VCF spec.
    return tuple(rec.id.split(";"))


def _pysam_filters(rec) -> tuple:
    return tuple(rec.filter.keys())


def _open_both(fixture_name: str):
    path = FIXTURE_DIR / fixture_name
    assert path.exists(), f"fixture missing: {path}"
    rf = rubam.VariantFile(str(path), "r")
    pf = pysam.VariantFile(str(path), "r")
    try:
        rrecs = list(rf)
        precs = list(pf)
    finally:
        rf.close()
        pf.close()
    return rrecs, precs


def _assert_core_tuple_match(rrecs, precs):
    assert len(rrecs) == len(precs), (
        f"record count mismatch: rubam={len(rrecs)} pysam={len(precs)}"
    )
    for i, (rr, pr) in enumerate(zip(rrecs, precs)):
        rt = (
            rr.reference_name,
            rr.position,
            rr.reference,
            tuple(rr.alternates or ()),
            tuple(rr.ids or ()),
            rr.quality,
            tuple(rr.filters or ()),
        )
        pt = (
            pr.chrom,
            pr.pos,
            pr.ref,
            tuple(pr.alts or ()),
            _pysam_ids(pr),
            pr.qual,
            _pysam_filters(pr),
        )
        assert rt == pt, f"record {i} core-tuple diff:\n  rubam={rt!r}\n  pysam={pt!r}"


def _assert_format_fields_match(rrecs, precs, fields=None):
    """Compare every FORMAT field for every sample, every record."""
    for i, (rr, pr) in enumerate(zip(rrecs, precs)):
        sample_names = tuple(pr.samples.keys())
        for sn in sample_names:
            ps = pr.samples[sn]
            rs = rr.samples[sn]
            pitems = dict(ps.items())
            keys = pitems.keys() if fields is None else fields
            for k in keys:
                pv = pitems.get(k)
                try:
                    rv = rs[k]
                except KeyError:
                    pytest.fail(
                        f"rec {i} sample {sn!r}: rubam missing FORMAT key {k!r}"
                    )
                assert rv == pv, (
                    f"rec {i} sample {sn!r} FORMAT {k}: rubam={rv!r} pysam={pv!r}"
                )


# ---------------------------------------------------------------------------
# Core-tuple parity: covers chrom/pos/ref/alts/ids/qual/filters.
# These are the fields the GIAB HG002 benchmark already validates; here we
# re-prove parity on every edge-case shape.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "fixture",
    [
        "snv_only.vcf",
        "mnv.vcf",
        "indels.vcf",
        "multi_allelic.vcf",
        "phased_gt.vcf",
        "missing_values.vcf",
        "multi_sample.vcf",
        "format_complex.vcf",
    ],
)
def test_core_tuple_parity(fixture):
    rrecs, precs = _open_both(fixture)
    _assert_core_tuple_match(rrecs, precs)


# ---------------------------------------------------------------------------
# Multi-allelic ALT tuple
# ---------------------------------------------------------------------------


def test_multi_allelic_alts_preserved():
    rrecs, precs = _open_both("multi_allelic.vcf")
    # First record: A -> G,T (2 alts); second record: C -> A,G,T (3 alts).
    assert rrecs[0].alternates == ("G", "T")
    assert rrecs[1].alternates == ("A", "G", "T")
    for rr, pr in zip(rrecs, precs):
        assert tuple(rr.alternates) == tuple(pr.alts)


# ---------------------------------------------------------------------------
# FORMAT fields: GT / DP / AD / GQ / PL across samples, including missing.
# These are the fields that pysam decodes per-sample; rubam exposes them via
# rec.samples[name][field].
# ---------------------------------------------------------------------------


def test_format_complex_all_fields():
    rrecs, precs = _open_both("format_complex.vcf")
    _assert_format_fields_match(rrecs, precs)


def test_multi_sample_format():
    rrecs, precs = _open_both("multi_sample.vcf")
    # Sanity: 5 samples declared.
    assert len(rrecs[0].samples) == 5
    _assert_format_fields_match(rrecs, precs)


def test_missing_values_format():
    """./. and '.' must decode identically in rubam vs pysam."""
    rrecs, precs = _open_both("missing_values.vcf")
    _assert_format_fields_match(rrecs, precs)


# ---------------------------------------------------------------------------
# Phased GT: rubam currently exposes the allele tuple but no `phased` flag.
# We assert GT tuple parity; the phased-flag check is xfail-marked.
# ---------------------------------------------------------------------------


def test_phased_gt_tuple_parity():
    """rubam returns the same (a, b) tuple as pysam, regardless of phasing."""
    rrecs, precs = _open_both("phased_gt.vcf")
    _assert_format_fields_match(rrecs, precs, fields=("GT", "DP"))


def test_phased_gt_flag_exposed():
    """**v0.3.2**: VariantSample.phased now exposed (closes the v3 xfail)."""
    rrecs, _ = _open_both("phased_gt.vcf")
    s_phased = rrecs[0].samples["SAMP1"]
    s_unphased = rrecs[0].samples["SAMP2"]
    assert hasattr(s_phased, "phased"), "rubam.VariantSample should expose .phased"
    assert s_phased.phased is True, "SAMP1 at pos=100 is 0|1 -> phased=True"
    assert s_unphased.phased is False, "SAMP2 at pos=100 is 0/1 -> phased=False"


# ---------------------------------------------------------------------------
# INFO field read-side access: rubam currently has only `set_info` (write),
# no `rec.info` getter. This is a known gap; xfail-marked, not skipped.
# ---------------------------------------------------------------------------


def test_info_read_side_parity():
    """**v0.3.2**: VariantRecord.info getter now exposed (closes v3 xfail).

    Compares the set of INFO keys with pysam. Value-level type parity is
    looser than the strict pysam dict equality the v3 xfail implied, because
    pysam returns numpy scalars for Number=1 numeric fields while rubam
    returns native Python int/float — we therefore compare keys + numerical
    near-equality.
    """
    rrecs, precs = _open_both("snv_only.vcf")
    for rr, pr in zip(rrecs, precs):
        rubam_info = dict(rr.info)
        pysam_info = dict(pr.info)
        assert set(rubam_info.keys()) == set(pysam_info.keys()), (
            f"INFO key mismatch: rubam={set(rubam_info)} pysam={set(pysam_info)}"
        )
        for key in rubam_info:
            r_val = rubam_info[key]
            p_val = pysam_info[key]
            if isinstance(r_val, (int, float)) and isinstance(p_val, (int, float)):
                assert abs(float(r_val) - float(p_val)) < 1e-3, f"{key}: {r_val} vs {p_val}"
            elif isinstance(r_val, tuple) and isinstance(p_val, (tuple, list)):
                assert len(r_val) == len(p_val)
                for a, b in zip(r_val, p_val):
                    if isinstance(a, (int, float)):
                        assert abs(float(a) - float(b)) < 1e-3
                    else:
                        assert a == b
            else:
                assert r_val == p_val, f"{key}: {r_val!r} vs {p_val!r}"


def test_info_number_A_round_trip():
    """**v0.3.2**: Number=A INFO fields parsed as tuple on the rubam side.

    Note: the underlying VCF stores AF as Type=Float which round-trips through
    f32; the tuple compares with tolerance rather than exact equality.
    """
    rrecs, _ = _open_both("multi_allelic.vcf")
    info = dict(rrecs[0].info)
    af = info.get("AF")
    assert isinstance(af, tuple) and len(af) == 2
    assert abs(af[0] - 0.3) < 1e-3
    assert abs(af[1] - 0.2) < 1e-3


# ---------------------------------------------------------------------------
# Out-of-scope per project policy: symbolic alleles, BND, polyploid.
# We do not author fixtures for them; documenting via xfail keeps the matrix
# auditable.
# ---------------------------------------------------------------------------


@pytest.mark.xfail(
    reason="Symbolic ALT (<DEL>, <INS>) parsing not yet validated. Will be "
    "added when rubam wires noodles symbolic-allele support.",
    strict=True,
    run=False,
)
def test_symbolic_alleles_not_yet_supported():
    raise NotImplementedError("symbolic alleles deferred")


@pytest.mark.xfail(
    reason="Breakend (BND) records not yet validated.",
    strict=True,
    run=False,
)
def test_breakends_not_yet_supported():
    raise NotImplementedError("BND deferred")


@pytest.mark.xfail(
    reason="Polyploid genotypes (>2 alleles per sample) not yet validated.",
    strict=True,
    run=False,
)
def test_polyploid_not_yet_supported():
    raise NotImplementedError("polyploid deferred")
