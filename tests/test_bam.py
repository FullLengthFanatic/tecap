"""Tests for BAM-level orchestration: accumulator invariants + thread equivalence."""


import pytest

from tecap.bam import analyse_bam, probe_strandedness
from tecap.constants import CAPTURED
from tecap.gtf import build_gene_index
from tecap.polya import build_polya_index


def _run(fx, threads, cb_tag=None, skip_antisense=False):
    gene_index, gene_records = build_gene_index(fx["gtf"])
    polya_index = build_polya_index(fx["polya"], slop=25,
                                    allowed_types={"TE", "DS", "DI", "EX"})
    return analyse_bam(
        bam_path      = fx["bam"],
        gene_index    = gene_index,
        gene_records  = gene_records,
        polya_index   = polya_index,
        cov_threshold = 0.5,
        min_mapq      = 0,
        threads       = threads,
        cb_tag        = cb_tag,
        skip_antisense = skip_antisense,
    )


def _probe(fx, **kw):
    gene_index, gene_records = build_gene_index(fx["gtf"])
    return probe_strandedness(fx["bam"], gene_index, gene_records, **kw)


def test_probe_strandedness_on_the_fixture(synthetic_fixture):
    # Every fixture read is sense to its gene (GA, GC forward; GB reverse).
    assert _probe(synthetic_fixture, min_mapq=0) == {
        "reads": len(synthetic_fixture["expected_buckets"]),
        "match_frac": 1.0, "stranded": True}


def test_probe_strandedness_on_half_flipped_reads(half_flipped_fixture):
    assert _probe(half_flipped_fixture, min_mapq=0) == {
        "reads": 12, "match_frac": 0.5, "stranded": False}


def test_probe_strandedness_stops_at_n(half_flipped_fixture):
    assert _probe(half_flipped_fixture, min_mapq=0, n=3)["reads"] == 3


def test_probe_strandedness_applies_the_mapq_filter(synthetic_fixture):
    # Every fixture read has MAPQ 60: no read qualifies, the input counts as
    # unstranded.
    assert _probe(synthetic_fixture, min_mapq=61) == {
        "reads": 0, "match_frac": None, "stranded": False}


@pytest.mark.parametrize("threads", [1, 2])
def test_antisense_read_is_skipped_on_stranded_input(antisense_fixture, threads):
    r = _run(antisense_fixture, threads, cb_tag="CB", skip_antisense=True)
    assert r["antisense_skipped"] == 1
    assert r["counts"][CAPTURED] == 3 and r["total"] == 12
    # Orientation counts every read with a gene, before the skip.
    assert r["orient"] == {"match": 12, "mismatch": 1}
    assert sum(r["site"].values()) == 12
    assert sum(r["utr_bin_total"]) == 12
    assert "CELL_A" not in r["per_cell"]
    assert r["per_gene"]["GA"][CAPTURED] == 2


@pytest.mark.parametrize("threads", [1, 2])
def test_antisense_read_is_classified_when_kept(antisense_fixture, threads):
    r = _run(antisense_fixture, threads, cb_tag="CB")
    assert r["antisense_skipped"] == 0
    assert r["counts"][CAPTURED] == 4 and r["total"] == 13
    assert r["orient"] == {"match": 12, "mismatch": 1}
    assert r["per_cell"]["CELL_A"] == {CAPTURED: 1}
    assert r["per_gene"]["GA"][CAPTURED] == 3


def test_accumulator_invariants(synthetic_fixture):
    r = _run(synthetic_fixture, threads=1)
    total = sum(r["counts"].values())
    assert r["total"] == total
    assert abs(sum(r["fractions"].values()) - 1.0) < 1e-9
    assert 0.0 <= (r["orient"]["match"] / (r["orient"]["match"] + r["orient"]["mismatch"])) <= 1.0


def test_single_vs_multi_thread_equivalence(synthetic_fixture):
    r1 = _run(synthetic_fixture, threads=1)
    r2 = _run(synthetic_fixture, threads=2)
    # Counts must match exactly
    assert r1["counts"] == r2["counts"]
    assert r1["pas_counts"] == r2["pas_counts"]
    assert r1["utr_bin_total"] == r2["utr_bin_total"]
    assert r1["utr_bin_captured"] == r2["utr_bin_captured"]
    assert r1["utr_bin_mecha_correct"] == r2["utr_bin_mecha_correct"]


def test_cb_tag_populated_when_present(synthetic_fixture):
    r = _run(synthetic_fixture, threads=1, cb_tag="CB")
    # GENE_A reads carry CB tags → per_cell should be non-empty.
    assert r["cb_tag_seen"] is True
    assert len(r["per_cell"]) >= 1


def test_cb_tag_absent_is_noop(synthetic_fixture):
    # Use a tag that's not set on any read.
    r = _run(synthetic_fixture, threads=1, cb_tag="XY")
    assert r["cb_tag_seen"] is False
    assert r["per_cell"] == {}


def test_per_gene_counts_present(synthetic_fixture):
    r = _run(synthetic_fixture, threads=1)
    # All three genes should have at least one read assigned.
    seen = set(r["per_gene"].keys())
    assert {"GA", "GB", "GC"}.issubset(seen)


def test_warn_if_unstranded_fires_below_threshold(caplog):
    import logging

    from tecap.bam import warn_if_unstranded
    log = logging.getLogger("tecap.test.orient")
    with caplog.at_level(logging.WARNING):
        warn_if_unstranded(0.49, log)
    assert "unstranded" in caplog.text.lower()


def test_warn_if_unstranded_silent_when_stranded(caplog):
    import logging

    from tecap.bam import warn_if_unstranded
    log = logging.getLogger("tecap.test.orient")
    with caplog.at_level(logging.WARNING):
        warn_if_unstranded(0.96, log)
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]


def test_warn_if_unstranded_silent_when_none(caplog):
    import logging

    from tecap.bam import warn_if_unstranded
    log = logging.getLogger("tecap.test.orient")
    with caplog.at_level(logging.WARNING):
        warn_if_unstranded(None, log)
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]
