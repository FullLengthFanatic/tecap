"""End-to-end CLI smoke tests via subprocess."""

import argparse
import json
import subprocess
import sys


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
