"""End-to-end CLI smoke tests via subprocess."""

import json
import os
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
