"""Tests for downstream base composition utilities."""

import pysam

from tecap.basecomp import analyse, downstream_window, pct_a, rc, summarise
from tecap.gtf import build_gene_index
from tecap.polya import build_polya_index


def test_rc():
    assert rc("ACGT")  == "ACGT"
    assert rc("AAAA")  == "TTTT"
    assert rc("NGCA")  == "TGCN"


def test_downstream_window_plus(synthetic_fixture):
    fasta = pysam.FastaFile(synthetic_fixture["fasta"])
    # We planted 20 A's at [9950,9970). A + strand read with 3' end = 9950
    # should see a window of 20 A's.
    seq = downstream_window(fasta, "chr_test", 9950, "+", 20)
    assert seq == "A" * 20
    assert pct_a(seq) == 1.0
    fasta.close()


def test_downstream_window_minus(synthetic_fixture):
    fasta = pysam.FastaFile(synthetic_fixture["fasta"])
    # On -, downstream_window returns RC of [pos-win, pos). At pos=9970, win=20
    # the slice is [9950,9970) = 20 A's, RC = 20 T's → 0% A.
    seq = downstream_window(fasta, "chr_test", 9970, "-", 20)
    assert seq == "T" * 20
    assert pct_a(seq) == 0.0
    fasta.close()


def test_downstream_window_out_of_bounds(synthetic_fixture):
    fasta = pysam.FastaFile(synthetic_fixture["fasta"])
    # Past contig end → short sequence → empty string returned
    seq = downstream_window(fasta, "chr_test", synthetic_fixture["contig_len"] - 5, "+", 20)
    assert seq == ""
    fasta.close()


def test_analyse_end_to_end(synthetic_fixture):
    fx = synthetic_fixture
    gene_index, gene_records = build_gene_index(fx["gtf"])
    polya_index = build_polya_index(fx["polya"], slop=25,
                                    allowed_types={"TE", "DS", "DI", "EX"})
    pct_a_by_bucket, bucket_counts = analyse(
        bam_path=fx["bam"], fasta_path=fx["fasta"],
        gene_index=gene_index, gene_records=gene_records,
        polya_index=polya_index,
        cov_threshold=0.5, min_mapq=0, window=20, threads=1,
    )
    summary = summarise(pct_a_by_bucket)
    assert "MechA_PAS+" in summary["buckets"]
    ma_pas = summary["buckets"]["MechA_PAS+"]
    # r_mechA_pas_1 has 3'end=9950, window 20 A's → pct_a=1.0.
    assert ma_pas["n"] >= 1
    assert ma_pas["median"] == 1.0
