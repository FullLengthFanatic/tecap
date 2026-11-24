"""Tests for BAM-level orchestration: accumulator invariants + thread equivalence."""


import pysam
import pytest

from tecap import bam as bam_mod
from tecap.bam import analyse_bam, choose_gene, gene_of_read, probe_strandedness
from tecap.constants import CAPTURED, MECH_B_ASPECI, MECH_B_EXON
from tecap.gtf import build_gene_index
from tecap.polya import build_polya_index
from tests.conftest import ANTISENSE_READ, CONTIG, CONTIG_LEN, REASSIGNED_READ


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


_HEADER = pysam.AlignmentHeader.from_dict(
    {"HD": {"VN": "1.6"}, "SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}]})


def _segment(start, cigar, reverse, seq=None):
    """A hand-built primary record on the fixture contig, MAPQ 60."""
    a = pysam.AlignedSegment(_HEADER)
    a.query_name = "r_hand"
    qlen = sum(n for op, n in cigar if op in (0, 1, 4, 7, 8))
    a.query_sequence = seq or "C" * qlen
    a.flag = 16 if reverse else 0
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = 60
    a.cigartuples = cigar
    return a


_ANTISENSE  = dict(start=ANTISENSE_READ[2], cigar=ANTISENSE_READ[3])
_REASSIGNED = dict(start=REASSIGNED_READ[2], cigar=REASSIGNED_READ[3],
                   seq=REASSIGNED_READ[4])


def test_gene_of_read_on_one_strand(reassign_fixture):
    gene_index, gene_records = build_gene_index(reassign_fixture["gtf"])
    genes = gene_index[CONTIG]
    read = _segment(reverse=True, **_REASSIGNED)
    assert gene_of_read(read, genes, gene_records, 0) == "GA"
    assert gene_of_read(read, genes, gene_records, 0, strand="+") == "GA"
    assert gene_of_read(read, genes, gene_records, 0, strand="-") == "GD"
    # Task 2.8's antisense read: no minus-strand gene overlaps it.
    read = _segment(reverse=True, **_ANTISENSE)
    assert gene_of_read(read, genes, gene_records, 0, strand="-") is None
    assert gene_of_read(read, genes, gene_records, 0, strand="+") == "GA"


@pytest.mark.parametrize("read, skip, want", [
    # Dropped by gene_of_read (no N in the CIGAR), with and without skip.
    (dict(start=3000, cigar=[(0, 600)], reverse=False), True,  (None, None)),
    (dict(start=3000, cigar=[(0, 600)], reverse=False), False, (None, None)),
    # Forward read on the + gene GA.
    (dict(reverse=False, **_ANTISENSE),  True,  ("GA", "sense")),
    (dict(reverse=False, **_ANTISENSE),  False, ("GA", "sense")),
    # Reverse on GA, unstranded: GA, as before Task 2.10.
    (dict(reverse=True,  **_ANTISENSE),  False, ("GA", "antisense")),
    (dict(reverse=True,  **_REASSIGNED), False, ("GA", "antisense")),
    # Reverse on GA, stranded: GD when a minus gene overlaps it, else skipped.
    (dict(reverse=True,  **_REASSIGNED), True,  ("GD", "reassigned")),
    (dict(reverse=True,  **_ANTISENSE),  True,  (None, "skipped")),
])
def test_choose_gene_outcomes(reassign_fixture, read, skip, want):
    gene_index, gene_records = build_gene_index(reassign_fixture["gtf"])
    got = choose_gene(_segment(**read), gene_index[CONTIG], gene_records, 0, skip)
    assert got == want


@pytest.mark.parametrize("threads", [1, 2])
def test_antisense_read_with_a_gene_on_its_strand_is_reassigned(reassign_fixture, threads):
    """Stranded: r_reassigned_D is classified against GD (MechB-exon, observed
    site), r_antisense_A has no gene on its strand and is skipped."""
    r = _run(reassign_fixture, threads, cb_tag="CB", skip_antisense=True)
    assert (r["antisense_skipped"], r["antisense_reassigned"]) == (1, 1)
    # Orientation is taken against the best-overlap gene on either strand.
    assert r["orient"] == {"match": 12, "mismatch": 2}
    assert r["total"] == 13
    assert (r["counts"][CAPTURED], r["counts"][MECH_B_EXON],
            r["counts"][MECH_B_ASPECI]) == (3, 2, 1)
    assert r["site"] == {"observed": 2, "fallback": 11}
    assert sum(r["utr_bin_total"]) == 13
    assert r["per_gene"]["GD"] == {MECH_B_EXON: 1}
    assert r["per_cell"]["CELL_D"] == {MECH_B_EXON: 1}
    assert "CELL_A" not in r["per_cell"]


@pytest.mark.parametrize("threads", [1, 2])
def test_antisense_reads_keep_their_gene_on_unstranded_input(reassign_fixture, threads):
    r = _run(reassign_fixture, threads, cb_tag="CB")
    assert (r["antisense_skipped"], r["antisense_reassigned"]) == (0, 0)
    assert r["orient"] == {"match": 12, "mismatch": 2}
    assert r["total"] == 14
    assert (r["counts"][CAPTURED], r["counts"][MECH_B_EXON],
            r["counts"][MECH_B_ASPECI]) == (4, 1, 2)
    assert r["site"] == {"observed": 1, "fallback": 13}
    assert "GD" not in r["per_gene"]
    assert r["per_cell"]["CELL_D"] == {MECH_B_ASPECI: 1}
    assert r["per_cell"]["CELL_A"] == {CAPTURED: 1}


@pytest.mark.parametrize("threads", [1, 2])
def test_antisense_read_is_skipped_on_stranded_input(antisense_fixture, threads):
    r = _run(antisense_fixture, threads, cb_tag="CB", skip_antisense=True)
    assert r["antisense_skipped"] == 1
    assert r["antisense_reassigned"] == 0
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


@pytest.mark.parametrize("anchor, offset", [("a_run_start", 0), ("alignment_end", 10)])
def test_worker_passes_the_site_to_classify_read(synthetic_fixture, monkeypatch, anchor, offset):
    fx = synthetic_fixture
    seen, real = {}, bam_mod.classify_read

    def spy(read, rec, polya_index, cov_threshold, read_3prime=None):
        seen[read.query_name] = read_3prime
        return real(read, rec, polya_index, cov_threshold, read_3prime=read_3prime)

    monkeypatch.setattr(bam_mod, "classify_read", spy)
    gene_index, gene_records = build_gene_index(fx["gtf"])
    polya_index = build_polya_index(fx["polya"], slop=25, allowed_types={"TE", "DS", "DI", "EX"})
    analyse_bam(bam_path=fx["bam"], gene_index=gene_index, gene_records=gene_records,
                polya_index=polya_index, cov_threshold=0.5, min_mapq=0, threads=1,
                site_anchor=anchor)
    assert seen["r_internal_kept"] == fx["internal_kept_site"] + offset


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
