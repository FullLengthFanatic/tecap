"""Tests for BAM-level orchestration: accumulator invariants + thread equivalence."""


from tecap.bam import analyse_bam
from tecap.gtf import build_gene_index
from tecap.polya import build_polya_index


def _run(fx, threads, cb_tag=None):
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
    )


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
