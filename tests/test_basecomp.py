"""Tests for downstream base composition utilities."""

import numpy as np
import pysam
import pytest

from tecap import basecomp
from tecap.bam import analyse_bam
from tecap.basecomp import analyse, assign_bucket, downstream_window, rc, summarise
from tecap.constants import DEFINITION_SOURCE, LABEL_WINDOW, LABELS
from tecap.gtf import build_gene_index
from tecap.polya import build_polya_index
from tecap.site import label_of
from tests.test_site import CONTIG_LEN, HEADER, L_VECTORS, build, expand


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
    assert seq.count("A") == 20
    fasta.close()


def test_downstream_window_minus(synthetic_fixture):
    fasta = pysam.FastaFile(synthetic_fixture["fasta"])
    # On -, downstream_window returns RC of [pos-win, pos). At pos=9970, win=20
    # the slice is [9950,9970) = 20 A's, RC = 20 T's → 0% A.
    seq = downstream_window(fasta, "chr_test", 9970, "-", 20)
    assert seq == "T" * 20
    assert seq.count("A") == 0
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
    a_count_by_bucket, bucket_counts, _, _ = analyse(
        bam_path=fx["bam"], fasta_path=fx["fasta"],
        gene_index=gene_index, gene_records=gene_records,
        polya_index=polya_index,
        cov_threshold=0.5, min_mapq=0, window=20, threads=1,
    )
    summary = summarise(a_count_by_bucket)
    assert "MechA_PAS+" in summary["buckets"]
    ma_pas = summary["buckets"]["MechA_PAS+"]
    # r_mechA_pas_1 has 3'end=9950, window 20 A's: 20 of 20, median 1.0.
    assert ma_pas["n"] >= 1
    assert ma_pas["median"] == 1.0


def test_analyse_reports_orientation(synthetic_fixture):
    fx = synthetic_fixture
    gene_index, gene_records = build_gene_index(fx["gtf"])
    polya_index = build_polya_index(fx["polya"], slop=25,
                                    allowed_types={"TE", "DS", "DI", "EX"})
    _, _, orient, _ = analyse(
        bam_path=fx["bam"], fasta_path=fx["fasta"],
        gene_index=gene_index, gene_records=gene_records,
        polya_index=polya_index,
        cov_threshold=0.5, min_mapq=0, window=20, threads=1,
    )
    # Every synthetic read is correctly oriented (GA forward/+, GB reverse/-).
    assert orient["match"] + orient["mismatch"] > 0
    assert orient["mismatch"] == 0


@pytest.mark.parametrize("anchor, offset, n_a", [
    ("a_run_start",   0, 10),    # window [P, P+20): the 10 planted A's, then C's
    ("alignment_end", 10, 0),    # window [P+10, P+30): C's only
])
def test_window_of_internal_kept_follows_the_anchor(synthetic_fixture, monkeypatch,
                                                    anchor, offset, n_a):
    fx = synthetic_fixture
    p = fx["internal_kept_site"]
    windows = {}

    def spy(fasta, chrom, pos, strand, win):
        windows[pos] = downstream_window(fasta, chrom, pos, strand, win)
        return windows[pos]

    monkeypatch.setattr(basecomp, "downstream_window", spy)
    gene_index, gene_records = build_gene_index(fx["gtf"])
    polya_index = build_polya_index(fx["polya"], slop=25,
                                    allowed_types={"TE", "DS", "DI", "EX"})
    analyse(
        bam_path=fx["bam"], fasta_path=fx["fasta"],
        gene_index=gene_index, gene_records=gene_records,
        polya_index=polya_index,
        cov_threshold=0.5, min_mapq=0, window=20, threads=1,
        site_anchor=anchor,
    )
    assert windows[p + offset].count("A") == n_a
    assert p + 10 - offset not in windows


@pytest.mark.parametrize("skip, skipped", [(True, 1), (False, 0)])
@pytest.mark.parametrize("threads", [1, 2])
def test_antisense_read_is_skipped_on_stranded_input(antisense_fixture, monkeypatch,
                                                     skip, skipped, threads):
    """r_antisense_A (Captured when kept, alignment end 9100) gets no window,
    no adenine count and no site when skip_antisense is set."""
    fx = antisense_fixture
    windows = []

    def spy(fasta, chrom, pos, strand, win):
        windows.append(pos)
        return downstream_window(fasta, chrom, pos, strand, win)

    monkeypatch.setattr(basecomp, "downstream_window", spy)
    gene_index, gene_records = build_gene_index(fx["gtf"])
    polya_index = build_polya_index(fx["polya"], slop=25,
                                    allowed_types={"TE", "DS", "DI", "EX"})
    a_count_by_bucket, bucket_counts, orient, site_stats = analyse(
        bam_path=fx["bam"], fasta_path=fx["fasta"],
        gene_index=gene_index, gene_records=gene_records,
        polya_index=polya_index,
        cov_threshold=0.5, min_mapq=0, window=20, threads=threads,
        skip_antisense=skip,
    )
    assert site_stats["antisense_skipped"] == skipped
    assert site_stats["antisense_reassigned"] == 0
    assert bucket_counts["Captured"] == 4 - skipped
    assert sorted(a_count_by_bucket["Captured"]) == [0] * (3 - skipped) + [10]
    assert orient == {"match": 12, "mismatch": 1}
    assert site_stats["observed"] + site_stats["fallback"] == 13 - skipped
    if threads == 1:                     # the spy is not seen by forked workers
        assert (9100 in windows) is not skip


@pytest.mark.parametrize("skip", [True, False])
@pytest.mark.parametrize("threads", [1, 2])
def test_reassigned_read_gets_the_site_of_its_own_gene(reassign_fixture, monkeypatch,
                                                       skip, threads):
    """Stranded: r_reassigned_D is bucketed against GD (MechB_exon) with its
    window at its observed site 2700 on minus; r_antisense_A is skipped.
    Unstranded: both stay on GA, r_reassigned_D with its window at 3600 on plus."""
    fx = reassign_fixture
    windows = []

    def spy(fasta, chrom, pos, strand, win):
        windows.append((pos, strand))
        return downstream_window(fasta, chrom, pos, strand, win)

    monkeypatch.setattr(basecomp, "downstream_window", spy)
    gene_index, gene_records = build_gene_index(fx["gtf"])
    polya_index = build_polya_index(fx["polya"], slop=25,
                                    allowed_types={"TE", "DS", "DI", "EX"})
    _, bucket_counts, orient, site_stats = analyse(
        bam_path=fx["bam"], fasta_path=fx["fasta"],
        gene_index=gene_index, gene_records=gene_records,
        polya_index=polya_index,
        cov_threshold=0.5, min_mapq=0, window=20, threads=threads,
        skip_antisense=skip,
    )
    assert orient == {"match": 12, "mismatch": 2}
    if skip:
        assert (site_stats["antisense_skipped"], site_stats["antisense_reassigned"]) == (1, 1)
        assert (bucket_counts["Captured"], bucket_counts["MechB_exon"],
                bucket_counts["MechB_aspecific"]) == (3, 2, 1)
        assert (site_stats["observed"], site_stats["fallback"]) == (2, 11)
        # The window [2680, 2700) on minus holds no adenine.
        assert site_stats["labels"]["MechB_exon"] == {
            "templated": 0, "nontemplated": 1, "ambiguous": 0}
    else:
        assert (site_stats["antisense_skipped"], site_stats["antisense_reassigned"]) == (0, 0)
        assert (bucket_counts["Captured"], bucket_counts["MechB_exon"],
                bucket_counts["MechB_aspecific"]) == (4, 1, 2)
        assert (site_stats["observed"], site_stats["fallback"]) == (1, 13)
        assert "MechB_exon" not in site_stats["labels"]
    if threads == 1:                     # the spy is not seen by forked workers
        assert ((2700, "-") in windows) is skip
        assert ((3600, "+") in windows) is not skip
        assert ((9100, "+") in windows) is not skip


@pytest.mark.parametrize("skip", [True, False])
def test_basecomp_and_classify_choose_the_same_gene(reassign_fixture, skip):
    """Same reads, same statuses and counts: orientation, skipped, reassigned,
    observed and fallback sites, and the bucket counts classify implies."""
    fx = reassign_fixture
    gene_index, gene_records = build_gene_index(fx["gtf"])
    polya_index = build_polya_index(fx["polya"], slop=25,
                                    allowed_types={"TE", "DS", "DI", "EX"})
    common = dict(gene_index=gene_index, gene_records=gene_records,
                  polya_index=polya_index, cov_threshold=0.5, min_mapq=0,
                  threads=1, skip_antisense=skip)
    r = analyse_bam(bam_path=fx["bam"], **common)
    _, bucket_counts, orient, site_stats = analyse(
        bam_path=fx["bam"], fasta_path=fx["fasta"], window=20, **common)

    assert orient == r["orient"]
    assert site_stats["antisense_skipped"] == r["antisense_skipped"]
    assert site_stats["antisense_reassigned"] == r["antisense_reassigned"]
    assert (site_stats["observed"], site_stats["fallback"]) == (
        r["site"]["observed"], r["site"]["fallback"])
    implied = {}
    for cat, n in r["counts"].items():
        split = r["pas_counts"].get(cat, {"pas": 0, "nopas": n})
        for has_pas, k in ((True, split["pas"]), (False, split["nopas"])):
            b = assign_bucket(cat, has_pas)
            if b is not None and k:
                implied[b] = implied.get(b, 0) + k
    assert bucket_counts == implied


def test_summarise_stamps_thresholds():
    out = summarise({})
    th = out["thresholds"]
    assert th["frac_ge_60"] == 0.60
    assert th["frac_30_50"] == [0.30, 0.50]
    assert th["frac_ge_40"] == 0.40
    assert th["hist_bins"] == 21


@pytest.mark.parametrize("anchor, want, a_counts", [
    # r_internal_kept is the only fixture read with an observed site; its
    # window [P, P+20) holds the 10 planted A's.
    ("a_run_start",   {"templated": 1, "nontemplated": 0, "ambiguous": 0}, [0, 0, 10]),
    # Every read is observed; the three Captured windows, r_internal_kept's
    # [P+10, P+30) among them, hold C's only.
    ("alignment_end", {"templated": 0, "nontemplated": 3, "ambiguous": 0}, [0, 0, 0]),
])
@pytest.mark.parametrize("threads", [1, 2])
def test_label_of_internal_kept_follows_the_anchor(synthetic_fixture, anchor, want, a_counts,
                                                   threads):
    fx = synthetic_fixture
    gene_index, gene_records = build_gene_index(fx["gtf"])
    polya_index = build_polya_index(fx["polya"], slop=25,
                                    allowed_types={"TE", "DS", "DI", "EX"})
    a_count_by_bucket, _, _, site_stats = analyse(
        bam_path=fx["bam"], fasta_path=fx["fasta"],
        gene_index=gene_index, gene_records=gene_records,
        polya_index=polya_index,
        cov_threshold=0.5, min_mapq=0, window=20, threads=threads,
        site_anchor=anchor,
    )
    assert site_stats["labels"]["Captured"] == want
    assert sorted(a_count_by_bucket["Captured"]) == a_counts
    assert all(type(k) is int for ks in a_count_by_bucket.values() for k in ks)


def _write_chrS(path, seq):
    with open(path, "w") as f:
        f.write(f">chrS\n{seq}\n")
    pysam.faidx(str(path))
    return str(path)


@pytest.mark.parametrize("name", [k for k, v in L_VECTORS.items() if v.site is not None])
@pytest.mark.parametrize("minus", [False, True])
def test_labels_of_the_shared_vectors(tmp_path, name, minus):
    """Appendix A contig of the row (its mirror on minus): the A count and the
    label of the window downstream of the table's site are the table's."""
    v = L_VECTORS[name]
    contig = "C" * CONTIG_LEN
    if v.planted:
        start, bases = v.planted[0], expand(v.planted[1])
        contig = contig[:start] + bases + contig[start + len(bases):]
    fasta = pysam.FastaFile(_write_chrS(tmp_path / "chrS.fa", rc(contig) if minus else contig))
    _, site = build(v, minus)
    window = downstream_window(fasta, "chrS", site, "-" if minus else "+", LABEL_WINDOW)
    fasta.close()
    if name == "L1_at_contig_end":
        assert window == ""
    else:
        assert window.count("A") == v.a_window
        assert label_of(window.count("A")) == v.label


def test_read_at_the_contig_end_gets_no_label(tmp_path, monkeypatch):    # Review Focus 3
    """L1_at_contig_end made spliced, and its mirror, both Captured: the run is
    found, the window falls off the contig, the read gets no adenine count and
    no label, and nothing crashes."""
    v = L_VECTORS["L1_at_contig_end"]._replace(pos=1790, cigar="10M100N90M10S")
    fasta = _write_chrS(tmp_path / "chrS.fa", "C" * CONTIG_LEN)
    gtf = tmp_path / "genes.gtf"
    gtf.write_text("".join(
        f'chrS\tt\texon\t{s + 1}\t{e}\t.\t{strand}\t.\tgene_id "{gid}";\n'
        for gid, strand, exons in (("GP", "+", [(1700, 1800), (1900, 2000)]),
                                   ("GM", "-", [(0, 100), (200, 300)]))
        for s, e in exons))
    bed = tmp_path / "polya.bed"
    bed.write_text("")
    bam_path = str(tmp_path / "end.bam")
    with pysam.AlignmentFile(bam_path, "wb", header=HEADER) as bam:
        for minus in (True, False):                  # coordinate order
            bam.write(build(v, minus)[0])
    pysam.index(bam_path)

    windows = []

    def spy(fasta, chrom, pos, strand, win):
        windows.append((pos, strand, downstream_window(fasta, chrom, pos, strand, win)))
        return windows[-1][2]

    monkeypatch.setattr(basecomp, "downstream_window", spy)
    gene_index, gene_records = build_gene_index(str(gtf))
    a_count_by_bucket, bucket_counts, _, site_stats = analyse(
        bam_path=bam_path, fasta_path=fasta,
        gene_index=gene_index, gene_records=gene_records,
        polya_index=build_polya_index(str(bed), slop=25),
        cov_threshold=0.5, min_mapq=0, window=LABEL_WINDOW, threads=1,
        site_anchor="a_run_start",
    )
    assert sorted(windows) == [(10, "-", ""), (CONTIG_LEN - 10, "+", "")]
    assert (site_stats["observed"], site_stats["fallback"]) == (2, 0)
    assert site_stats["labels"] == {}
    assert dict(a_count_by_bucket) == {} and bucket_counts == {}


def test_old_statistics_keep_their_values():
    counts = {"Captured": [0, 3, 6, 7, 8, 9, 12, 20]}
    old = np.asarray([k / 20 for k in counts["Captured"]], dtype=np.float32)
    b = summarise(counts)["buckets"]["Captured"]
    assert b["median"] == float(np.median(old)) and b["mean"] == float(np.mean(old))
    assert b["frac_ge_60"] == float((old >= 0.60).mean()) and b["frac_ge_40"] == float((old >= 0.40).mean())
    assert b["frac_30_50"] == float(((old >= 0.30) & (old < 0.50)).mean())
    assert b["hist"] == np.histogram(old, bins=np.linspace(0.0, 1.0, 21))[0].tolist()


def test_label_fractions_come_from_integer_counts():
    labels = {"Captured": {"templated": 3, "nontemplated": 2, "ambiguous": 1}}
    b = summarise({"Captured": [0, 3, 6, 8, 9, 20, 20]}, labels=labels)["buckets"]["Captured"]
    assert (b["n"], b["n_site"]) == (7, 6)
    assert (b["n_templated"], b["n_nontemplated"], b["n_ambiguous"]) == (3, 2, 1)
    assert b["frac_templated"] + b["frac_nontemplated"] + b["frac_ambiguous"] == pytest.approx(1.0)
    assert b["frac_templated"] == 0.5


def test_labels_are_null_when_the_window_is_not_20():        # Review Focus 4
    b = summarise({"Captured": [0, 15, 30]}, window=30, labels={"Captured": {"templated": 1, "nontemplated": 1, "ambiguous": 1}})["buckets"]["Captured"]
    assert b["frac_templated"] is None and b["median"] == 0.5


def test_thresholds_carry_the_definition():
    th = summarise({})["thresholds"]
    assert (th["templated_min"], th["nontemplated_max"], th["label_window"]) == (8, 4, 20)
    assert (th["site_run_min"], th["site_run_tolerance"]) == (10, 1)
    assert th["frac_ge_60"] == 0.60 and th["hist_bins"] == 21


def test_label_fields_are_null_without_labels_and_zero_without_sites():
    b = summarise({"Captured": [0, 10]})["buckets"]["Captured"]
    assert b["n_site"] is None
    assert all(b[f"{p}_{x}"] is None for p in ("n", "frac") for x in LABELS)
    b = summarise({"Captured": [0, 10]}, labels={})["buckets"]["Captured"]
    assert (b["n_site"], b["n_templated"], b["frac_templated"]) == (0, 0, None)


def test_summary_names_the_definition_source():
    assert summarise({})["definition_source"] == DEFINITION_SOURCE
