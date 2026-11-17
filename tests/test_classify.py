"""Tests for per-read classification: every read should land in its bucket."""

import pysam

from tecap.classify import classify_read
from tecap.constants import (
    INTERNAL_PRIME_TE_CDS,
    MECH_A_CORRECT,
    MECH_B_EXON,
    MECH_C,
)
from tecap.gtf import build_gene_index
from tecap.polya import build_polya_index


def _load_indices(fx):
    gene_index, gene_records = build_gene_index(fx["gtf"])
    polya_index = build_polya_index(fx["polya"], slop=25,
                                    allowed_types={"TE", "DS", "DI", "EX"})
    return gene_index, gene_records, polya_index


def _find_gene(read, gene_index, gene_records):
    chrom = read.reference_name
    hits  = gene_index[chrom][read.reference_start:read.reference_end]
    best_gid, best_ovlp = None, 0
    for iv in hits:
        gid = iv.data
        ovlp = (min(read.reference_end, iv.end) -
                max(read.reference_start, iv.begin))
        if ovlp > best_ovlp:
            best_ovlp, best_gid = ovlp, gid
    return best_gid


def test_every_read_classifies_as_expected(synthetic_fixture):
    gene_index, gene_records, polya_index = _load_indices(synthetic_fixture)

    bam = pysam.AlignmentFile(synthetic_fixture["bam"], "rb")
    expected = synthetic_fixture["expected_buckets"]
    has_pas_expect = synthetic_fixture["expected_has_pas"]

    seen = {}
    for read in bam.fetch(until_eof=True):
        if (read.is_unmapped or read.is_secondary or read.is_supplementary
                or not read.cigartuples):
            continue
        if not any(op == 3 for op, _ in read.cigartuples):
            continue
        gid = _find_gene(read, gene_index, gene_records)
        assert gid is not None, f"read {read.query_name} matched no gene"
        rec = gene_records[gid]
        cat, cov, _, has_pas = classify_read(read, rec, polya_index, 0.5)
        seen[read.query_name] = (cat, has_pas)

        assert cat == expected[read.query_name], (
            f"read {read.query_name}: got {cat!r}, expected {expected[read.query_name]!r}"
        )
        if read.query_name in has_pas_expect:
            assert has_pas == has_pas_expect[read.query_name], (
                f"read {read.query_name}: has_pas={has_pas}, "
                f"expected {has_pas_expect[read.query_name]}"
            )
    bam.close()

    assert set(seen) == set(expected), (
        f"missing reads: expected {set(expected)}, got {set(seen)}"
    )


def _read_named(bam, name):
    return next(r for r in bam.fetch(until_eof=True) if r.query_name == name)


def test_read_3prime_replaces_the_alignment_end(synthetic_fixture):
    _, gene_records, polya_index = _load_indices(synthetic_fixture)
    p = synthetic_fixture["internal_kept_site"]
    with pysam.AlignmentFile(synthetic_fixture["bam"], "rb") as bam:
        read = _read_named(bam, "r_internal_kept")
        rec = gene_records["GA"]
        _, cov_end, r3p_end, _ = classify_read(read, rec, polya_index, 0.5)
        _, cov_site, r3p_site, _ = classify_read(read, rec, polya_index, 0.5,
                                                 read_3prime=p)
    assert r3p_end == p + 10
    assert r3p_site == p
    assert cov_site == cov_end


def test_read_3prime_drives_every_position_test(synthetic_fixture):
    # r_mechA_intutr_1: alignment end 8300 (MechA-internalUTR), cov 0.40.
    _, gene_records, polya_index = _load_indices(synthetic_fixture)
    with pysam.AlignmentFile(synthetic_fixture["bam"], "rb") as bam:
        read = _read_named(bam, "r_mechA_intutr_1")
        rec = gene_records["GA"]
        got = {r3p: classify_read(read, rec, polya_index, 0.5, read_3prime=r3p)
               for r3p in (9950, 7500, 3200, 10500)}
    assert got[9950][0] == MECH_A_CORRECT and got[9950][3] is True
    assert got[7500][0] == INTERNAL_PRIME_TE_CDS
    assert got[3200][0] == MECH_B_EXON
    assert got[10500][0] == MECH_C
    assert {r3p: g[2] for r3p, g in got.items()} == {r: r for r in got}
