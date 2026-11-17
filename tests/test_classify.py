"""Tests for per-read classification: every read should land in its bucket."""

import pysam
import pytest
from intervaltree import IntervalTree

from tecap.classify import classify_read
from tecap.constants import (
    INTERNAL_PRIME_TE_CDS,
    MECH_A_CORRECT,
    MECH_B_APA,
    MECH_B_ASPECI,
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


# Exon and polyA-cluster edges. read_3prime on plus genes is the first base
# after the read; on minus genes it is the last cDNA base. The minus gene
# mirrors the plus gene: base p -> 11999 - p, interval [s, e) -> [12000 - e, 12000 - s).
EDGE_CHROM = "chr_edge"

# Plus: upstream exons [1000,1500) and [3000,3500), TE [7000,10000), CDS ends at 8000.
EDGE_PLUS = {
    "chrom": EDGE_CHROM, "strand": "+",
    "te_start": 7000, "te_end": 10000, "te_length": 3000,
    "cds_end_genomic": 8000,
    "upstream_exon_tree": IntervalTree.from_tuples([(1000, 1500, 1), (3000, 3500, 1)]),
}
# Minus: upstream exons [10500,11000) and [8500,9000), TE [2000,5000), CDS ends at 4000.
EDGE_MINUS = {
    "chrom": EDGE_CHROM, "strand": "-",
    "te_start": 2000, "te_end": 5000, "te_length": 3000,
    "cds_end_genomic": 4000,
    "upstream_exon_tree": IntervalTree.from_tuples([(8500, 9000, 1), (10500, 11000, 1)]),
}

# PAS+ clusters in the intron upstream of the TE; with slop 25 the plus one
# covers [4975,5125) (last base 5124) and the minus one [6875,7025) (last
# base on the transcript strand 6875). The 3' UTR clusters cover [8975,9125)
# on plus (last base 9124) and [2875,3025) on minus (last base 2875).
EDGE_CLUSTERS = [
    (EDGE_CHROM, 5000, 5100, "plus_up",   "0", "+", "0", "0", "0", "DS", "AATAAA"),
    (EDGE_CHROM, 6900, 7000, "minus_up",  "0", "-", "0", "0", "0", "DS", "AATAAA"),
    (EDGE_CHROM, 9000, 9100, "plus_utr",  "0", "+", "0", "0", "0", "TE", "AATAAA"),
    (EDGE_CHROM, 2900, 3000, "minus_utr", "0", "-", "0", "0", "0", "TE", "AATAAA"),
]


@pytest.fixture
def edge_polya_index(tmp_path):
    bed = tmp_path / "edge_polya.bed"
    bed.write_text("".join("\t".join(map(str, row)) + "\n" for row in EDGE_CLUSTERS))
    return build_polya_index(str(bed), slop=25, allowed_types={"TE", "DS", "DI", "EX"})


def _edge_read(start, end, strand):
    read = pysam.AlignedSegment()
    read.query_name = f"edge_{strand}_{start}_{end}"
    read.flag = 16 if strand == "-" else 0
    read.reference_start = start
    read.cigartuples = [(0, end - start)]
    return read


# (strand, read start, read end, last cDNA base, expected bucket)
EDGE_CASES = [
    # last cDNA base is the last base of the upstream exon
    ("+", 3000, 3500, 3499, MECH_B_EXON),
    ("-", 8500, 9000, 8500, MECH_B_EXON),
    # last cDNA base is the first intronic base after that exon
    ("+", 3000, 3501, 3500, MECH_B_ASPECI),
    ("-", 8499, 9000, 8499, MECH_B_ASPECI),
    # last cDNA base is the last position of the polyA cluster after slop
    ("+", 4900, 5125, 5124, MECH_B_APA),
    ("-", 6875, 7200, 6875, MECH_B_APA),
    # one base further: not APA
    ("+", 4900, 5126, 5125, MECH_B_ASPECI),
    ("-", 6874, 7200, 6874, MECH_B_ASPECI),
    # leading edges: last cDNA base one base before the exon, and one base
    # before the polyA cluster interval after slop
    ("+", 2500, 3000, 2999, MECH_B_ASPECI),
    ("-", 9000, 9500, 9000, MECH_B_ASPECI),
    ("+", 4900, 4975, 4974, MECH_B_ASPECI),
    ("-", 7025, 7200, 7025, MECH_B_ASPECI),
    # 3' UTR branch: last cDNA base is the last position of a UTR cluster
    # after slop (coverage 625 / 3000, below the 0.5 threshold)
    ("+", 8500, 9125, 9124, MECH_A_CORRECT),
    ("-", 2875, 3500, 2875, MECH_A_CORRECT),
]


@pytest.mark.parametrize("strand,start,end,last_base,expected", EDGE_CASES)
def test_position_tests_use_the_last_cdna_base(edge_polya_index, strand, start, end,
                                               last_base, expected):
    rec = EDGE_PLUS if strand == "+" else EDGE_MINUS
    read = _edge_read(start, end, strand)
    boundary = end if strand == "+" else start
    assert boundary == (last_base + 1 if strand == "+" else last_base)

    by_alignment = classify_read(read, rec, edge_polya_index, 0.5)
    by_site = classify_read(read, rec, edge_polya_index, 0.5, read_3prime=boundary)

    assert by_alignment[0] == expected
    assert by_alignment[2] == boundary
    assert by_alignment[3] is (expected in (MECH_B_APA, MECH_A_CORRECT))
    assert by_site == by_alignment
