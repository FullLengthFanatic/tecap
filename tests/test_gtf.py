"""Tests for GTF parsing and gene-record construction."""

from tecap.gtf import build_gene_index


def test_gene_coordinates(synthetic_fixture):
    gene_index, gene_records = build_gene_index(synthetic_fixture["gtf"])

    assert set(gene_records) == {"GA", "GB", "GC"}

    # GENE_A: + strand, TE = [7000, 10000], CDS ends at 8000, UTR = 2000
    ga = gene_records["GA"]
    assert ga["strand"] == "+"
    assert ga["te_start"] == 7000
    assert ga["te_end"]   == 10000
    assert ga["te_length"] == 3000
    assert ga["cds_end_genomic"] == 8000
    assert ga["utr_length"] == 2000

    # GENE_B: - strand, TE = [12000, 14000], cds_end_genomic = 12500, UTR = 500
    gb = gene_records["GB"]
    assert gb["strand"] == "-"
    assert gb["te_start"] == 12000
    assert gb["te_end"]   == 14000
    assert gb["te_length"] == 2000
    assert gb["cds_end_genomic"] == 12500
    assert gb["utr_length"] == 500

    # GENE_C: non-coding
    gc = gene_records["GC"]
    assert gc["strand"] == "+"
    assert gc["te_start"] == 22000
    assert gc["te_end"]   == 24000
    assert gc["cds_end_genomic"] is None
    assert gc["utr_length"] == 2000


def test_gene_index_overlap(synthetic_fixture):
    gene_index, _ = build_gene_index(synthetic_fixture["gtf"])
    hits = gene_index["chr_test"][5000:5001]
    assert len(hits) == 1
    assert next(iter(hits)).data == "GA"

    hits = gene_index["chr_test"][16000:16001]
    assert len(hits) == 1
    assert next(iter(hits)).data == "GB"

    hits = gene_index["chr_test"][21000:21001]
    assert len(hits) == 1
    assert next(iter(hits)).data == "GC"


def test_upstream_exon_tree_excludes_te(synthetic_fixture):
    _, gene_records = build_gene_index(synthetic_fixture["gtf"])
    ga_upstream = gene_records["GA"]["upstream_exon_tree"]
    # TE (7000-10000) should NOT be in upstream tree.
    assert len(ga_upstream[7500:7501]) == 0
    # Upstream exons ARE in it.
    assert len(ga_upstream[1200:1201]) == 1
    assert len(ga_upstream[3200:3201]) == 1
