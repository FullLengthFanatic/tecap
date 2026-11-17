"""End-to-end CLI smoke tests via subprocess."""

import argparse
import json
import subprocess
import sys

import pysam
import pytest

from tecap.constants import BUCKETS, DEFINITION_SOURCE, SITE_ANCHOR_DEFAULT


def _run_cli(args, cwd=None):
    cmd = [sys.executable, "-m", "tecap", *args]
    res = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    return res


def test_version_flag():
    res = _run_cli(["--version"])
    assert res.returncode == 0
    assert "tecap" in res.stdout


def test_drna_refused():
    res = _run_cli(["classify", "--bam", "x", "--gtf", "x",
                    "--polya-sites", "x", "--platform", "drna-ont"])
    assert res.returncode != 0
    assert "Direct-RNA" in res.stderr or "drna" in res.stderr.lower()


def test_classify_end_to_end(synthetic_fixture, tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    res = _run_cli([
        "classify",
        "--bam", synthetic_fixture["bam"],
        "--gtf", synthetic_fixture["gtf"],
        "--polya-sites", synthetic_fixture["polya"],
        "--sample", "TST",
        "--out-dir", str(out),
        "--platform", "cdna-pacbio",
        "--per-gene-table",
        "--cell-barcode-tag", "CB",
    ])
    assert res.returncode == 0, f"stderr: {res.stderr}"

    j = out / "TST_terminal_exon.json"
    assert j.exists()
    data = json.loads(j.read_text())
    assert data["schema_version"] == "0.1"
    assert data["platform"] == "cdna-pacbio"
    assert data["sample"] == "TST"
    assert data["total"] > 0

    # Per-gene TSV
    tsv = out / "TST_per_gene.tsv"
    assert tsv.exists()
    lines = tsv.read_text().strip().splitlines()
    assert lines[0].startswith("gene_id\tchrom\tstrand\tte_length\tutr_length\t")
    assert len(lines) >= 2

    # Per-cell metrics populated because GA reads carry CB.
    assert "per_cell_metrics" in data
    assert len(data["per_cell_metrics"]) >= 1


def test_basecomp_end_to_end(synthetic_fixture, tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    res = _run_cli([
        "basecomp",
        "--bam", synthetic_fixture["bam"],
        "--gtf", synthetic_fixture["gtf"],
        "--polya-sites", synthetic_fixture["polya"],
        "--fasta", synthetic_fixture["fasta"],
        "--sample", "TST",
        "--out-dir", str(out),
        "--platform", "cdna-pacbio",
    ])
    assert res.returncode == 0, f"stderr: {res.stderr}"

    j = out / "TST_basecomp.json"
    assert j.exists()
    data = json.loads(j.read_text())
    assert data["schema_version"] == "0.1"
    assert data["sample"] == "TST"
    assert data["window"] == 20


def test_basecomp_json_is_self_describing(synthetic_fixture, tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    res = _run_cli([
        "basecomp",
        "--bam", synthetic_fixture["bam"],
        "--gtf", synthetic_fixture["gtf"],
        "--polya-sites", synthetic_fixture["polya"],
        "--fasta", synthetic_fixture["fasta"],
        "--sample", "TST",
        "--out-dir", str(out),
        "--platform", "cdna-pacbio",
    ])
    assert res.returncode == 0, f"stderr: {res.stderr}"

    data = json.loads((out / "TST_basecomp.json").read_text())
    # Orientation is recorded; synthetic reads are all correctly oriented.
    assert data["orient_match_frac"] == 1.0
    # The regime cutoffs that produced frac_* are stamped into the file.
    assert data["thresholds"]["frac_ge_60"] == 0.60
    assert data["thresholds"]["frac_30_50"] == [0.30, 0.50]


def test_basecomp_json_carries_the_labels(synthetic_fixture, tmp_path):
    fx = synthetic_fixture
    res = _run_cli(["basecomp", "--bam", fx["bam"], "--gtf", fx["gtf"],
                    "--polya-sites", fx["polya"], "--fasta", fx["fasta"],
                    "--sample", "TST", "--out-dir", str(tmp_path)])
    assert res.returncode == 0, f"stderr: {res.stderr}"
    assert "label fields are null" not in res.stderr
    data = json.loads((tmp_path / "TST_basecomp.json").read_text())
    assert data["definition_source"] == DEFINITION_SOURCE
    assert data["thresholds"]["templated_min"] == 8
    # Default anchor: r_internal_kept is the one read with an observed site.
    cap = data["buckets"]["Captured"]
    assert (cap["n"], cap["n_site"], cap["n_templated"], cap["frac_templated"]) == (3, 1, 1, 1.0)
    assert data["buckets"]["MechA_PAS+"]["n_site"] == 0
    assert data["site_counts"] == {"observed": 1, "fallback": len(fx["expected_buckets"]) - 1}


def test_basecomp_window_other_than_20(synthetic_fixture, tmp_path):     # Review Focus 4
    """Bands computed on the 30 nt window as before; label fields null; warned."""
    fx = synthetic_fixture
    res = _run_cli(["basecomp", "--bam", fx["bam"], "--gtf", fx["gtf"],
                    "--polya-sites", fx["polya"], "--fasta", fx["fasta"],
                    "--sample", "TST", "--out-dir", str(tmp_path), "--window", "30"])
    assert res.returncode == 0, f"stderr: {res.stderr}"
    assert "WARNING" in res.stderr and "label fields are null" in res.stderr
    data = json.loads((tmp_path / "TST_basecomp.json").read_text())
    assert data["window"] == 30
    # r_mechA_pas_1: window [9950, 9980) holds 20 planted A's then C's.
    ma_pas = data["buckets"]["MechA_PAS+"]
    assert ma_pas["median"] == pytest.approx(20 / 30) and ma_pas["frac_ge_60"] == 1.0
    for b in BUCKETS:
        d = data["buckets"][b]
        assert all(d[k] is None for k in ("n_site", "n_templated", "n_nontemplated",
                                          "n_ambiguous", "frac_templated",
                                          "frac_nontemplated", "frac_ambiguous")), b


@pytest.mark.parametrize("anchor", [None, "alignment_end", "a_run_start"])
def test_site_anchor_and_counts_in_both_jsons(synthetic_fixture, tmp_path, anchor):
    fx = synthetic_fixture
    n_reads = len(fx["expected_buckets"])
    used = anchor or SITE_ANCHOR_DEFAULT
    # Only r_internal_kept carries a terminal A run; the other reads are C's.
    want = ({"observed": n_reads, "fallback": 0} if used == "alignment_end"
            else {"observed": 1, "fallback": n_reads - 1})
    common = ["--bam", fx["bam"], "--gtf", fx["gtf"], "--polya-sites", fx["polya"],
              "--sample", "TST", "--out-dir", str(tmp_path)]
    if anchor is not None:
        common += ["--site-anchor", anchor]

    for cmd, extra, json_name in (
        ("classify", [], "TST_terminal_exon.json"),
        ("basecomp", ["--fasta", fx["fasta"]], "TST_basecomp.json"),
    ):
        res = _run_cli([cmd, *common, *extra])
        assert res.returncode == 0, f"{cmd} stderr: {res.stderr}"
        data = json.loads((tmp_path / json_name).read_text())
        assert data["site_anchor"] == used, cmd
        assert data["site_counts"] == want, cmd
        if cmd == "classify":
            assert data["total"] == n_reads


@pytest.mark.parametrize("threads", ["1", "2"])
def test_records_without_sequence_fall_back(synthetic_fixture, tmp_path, threads):
    """Primary spliced records stored without SEQ: both subcommands finish under
    the default anchor, classify them as before and count them as fallback."""
    fx = synthetic_fixture
    bam_path = str(tmp_path / "noseq.bam")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": fx["contig"], "LN": fx["contig_len"]}]}
    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for name, gid, start, cigar, *_ in fx["reads"]:
            if name not in ("r_captured_1", "r_captured_B"):   # one per strand
                continue
            a = pysam.AlignedSegment()
            a.query_name      = name
            a.flag            = 16 if gid == "GB" else 0
            a.reference_id    = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigartuples     = cigar
            bam.write(a)
    pysam.index(bam_path)

    common = ["--bam", bam_path, "--gtf", fx["gtf"], "--polya-sites", fx["polya"],
              "--sample", "NOSEQ", "--out-dir", str(tmp_path), "--threads", threads]
    res = _run_cli(["classify", *common])
    assert res.returncode == 0, f"classify stderr: {res.stderr}"
    data = json.loads((tmp_path / "NOSEQ_terminal_exon.json").read_text())
    assert data["site_anchor"] == SITE_ANCHOR_DEFAULT == "a_run_start"
    assert data["site_counts"] == {"observed": 0, "fallback": 2}
    assert data["counts"]["Captured (>=50%)"] == data["total"] == 2

    res = _run_cli(["basecomp", *common, "--fasta", fx["fasta"]])
    assert res.returncode == 0, f"basecomp stderr: {res.stderr}"
    data = json.loads((tmp_path / "NOSEQ_basecomp.json").read_text())
    assert data["site_counts"] == {"observed": 0, "fallback": 2}
    assert data["buckets"]["Captured"]["n"] == 2


def test_unknown_site_anchor_is_refused(synthetic_fixture, tmp_path):
    res = _run_cli(["classify", "--bam", synthetic_fixture["bam"],
                    "--gtf", synthetic_fixture["gtf"],
                    "--polya-sites", synthetic_fixture["polya"],
                    "--out-dir", str(tmp_path), "--site-anchor", "alignment-end"])
    assert res.returncode == 2
    assert "--site-anchor: invalid choice" in res.stderr


def test_explain_lists_all_mechanisms():
    res = _run_cli(["explain"])
    assert res.returncode == 0, res.stderr
    assert "Captured" in res.stdout
    assert "MechA-correct" in res.stdout
    assert "MechC" in res.stdout


def test_explain_single_mechanism_json():
    res = _run_cli(["explain", "--mechanism", "MechA-correct", "--format", "json"])
    assert res.returncode == 0, res.stderr
    payload = json.loads(res.stdout)
    assert "classify" in payload
    keys = list(payload["classify"].keys())
    assert any("MechA-correct" in k for k in keys)


def test_explain_unknown_mechanism_errors():
    res = _run_cli(["explain", "--mechanism", "not-a-mechanism"])
    assert res.returncode != 0


def test_report_subcommand(synthetic_fixture, tmp_path):
    """End-to-end: classify -> report HTML."""
    out = tmp_path / "out"
    out.mkdir()
    res = _run_cli([
        "classify",
        "--bam", synthetic_fixture["bam"],
        "--gtf", synthetic_fixture["gtf"],
        "--polya-sites", synthetic_fixture["polya"],
        "--sample", "TST",
        "--out-dir", str(out),
        "--platform", "cdna-pacbio",
    ])
    assert res.returncode == 0, res.stderr

    html_path = tmp_path / "report.html"
    res = _run_cli([
        "report",
        "--classify-json", str(out / "TST_terminal_exon.json"),
        "--out-html", str(html_path),
    ])
    assert res.returncode == 0, res.stderr
    text = html_path.read_text()
    assert "<html" in text
    assert "Captured" in text
    assert "Mechanism legend" in text


def test_report_multi_sample(synthetic_fixture, tmp_path):
    """Multi-sample report accepts space-separated --classify-json paths."""
    out_a = tmp_path / "out_a"
    out_b = tmp_path / "out_b"
    out_a.mkdir()
    out_b.mkdir()
    for sample, out in (("A", out_a), ("B", out_b)):
        res = _run_cli([
            "classify",
            "--bam", synthetic_fixture["bam"],
            "--gtf", synthetic_fixture["gtf"],
            "--polya-sites", synthetic_fixture["polya"],
            "--sample", sample,
            "--out-dir", str(out),
            "--platform", "cdna-pacbio",
        ])
        assert res.returncode == 0, res.stderr

    html_path = tmp_path / "compare.html"
    res = _run_cli([
        "report",
        "--classify-json",
        str(out_a / "A_terminal_exon.json"),
        str(out_b / "B_terminal_exon.json"),
        "--out-html", str(html_path),
    ])
    assert res.returncode == 0, res.stderr
    text = html_path.read_text()
    assert "<html" in text


def test_genome_triggers_auto_download(monkeypatch, synthetic_fixture, tmp_path):
    """--genome with no explicit paths should call fetch_polya_atlas /
    fetch_gencode_gtf instead of erroring."""
    from tecap import cli, download

    cache = tmp_path / "cache"

    def fake_fetch_polya(genome, out_dir, **kw):
        return synthetic_fixture["polya"], "deadbeef"

    def fake_fetch_gtf(genome, version, out_dir, **kw):
        return synthetic_fixture["gtf"], "cafef00d"

    monkeypatch.setattr(download, "fetch_polya_atlas", fake_fetch_polya)
    monkeypatch.setattr(download, "fetch_gencode_gtf", fake_fetch_gtf)

    args = argparse.Namespace(
        polya_sites=None, gtf=None,
        genome="GRCh38", gtf_version=45,
        ref_cache=str(cache),
    )
    cli._resolve_references(args)
    assert args.polya_sites == synthetic_fixture["polya"]
    assert args.gtf == synthetic_fixture["gtf"]


def test_genome_grcm39_polya_is_rejected(monkeypatch, tmp_path):
    """GRCm39 has no PolyASite atlas; without --polya-sites this must error."""
    from tecap import cli

    args = argparse.Namespace(
        polya_sites=None, gtf="x.gtf",
        genome="GRCm39", gtf_version=35,
        ref_cache=str(tmp_path),
    )
    try:
        cli._resolve_references(args)
    except SystemExit as e:
        assert e.code == 2
    else:
        raise AssertionError("expected SystemExit")
