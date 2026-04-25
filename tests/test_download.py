"""Tests for the download helper. Network is mocked; we check argument
plumbing, SHA256 verification, and the resume-skip behaviour."""

import hashlib

import pytest

from tecap import download


def _fake_download(content):
    def _impl(url, dst):
        with open(dst, "wb") as f:
            f.write(content)
        return dst
    return _impl


def test_fetch_polya_atlas_writes_and_returns_sha(monkeypatch, tmp_path):
    payload = b"fake polyA cluster bed"
    monkeypatch.setattr(download, "_download", _fake_download(payload))
    path, sha = download.fetch_polya_atlas("GRCh38", str(tmp_path))
    assert path.endswith(".bed.gz")
    assert sha == hashlib.sha256(payload).hexdigest()


def test_fetch_polya_atlas_unknown_genome_raises(tmp_path):
    with pytest.raises(ValueError, match="No built-in URL"):
        download.fetch_polya_atlas("GRCm38", str(tmp_path))


def test_fetch_polya_atlas_sha_mismatch(monkeypatch, tmp_path):
    monkeypatch.setattr(download, "_download", _fake_download(b"abc"))
    with pytest.raises(ValueError, match="SHA256 mismatch"):
        download.fetch_polya_atlas(
            "GRCh38", str(tmp_path), expected_sha256="0" * 64,
        )


def test_fetch_polya_atlas_skips_existing(monkeypatch, tmp_path):
    payload = b"existing"
    target = tmp_path / "atlas.clusters.3.0.GRCh38.bed.gz"
    target.write_bytes(payload)

    called = {"n": 0}

    def _should_not_run(url, dst):
        called["n"] += 1
        return dst

    monkeypatch.setattr(download, "_download", _should_not_run)
    path, sha = download.fetch_polya_atlas("GRCh38", str(tmp_path))
    assert called["n"] == 0
    assert sha == hashlib.sha256(payload).hexdigest()


def test_fetch_gencode_gtf_url_for_human(monkeypatch, tmp_path):
    seen = {}

    def _capture(url, dst):
        seen["url"] = url
        with open(dst, "wb") as f:
            f.write(b"x")
        return dst

    monkeypatch.setattr(download, "_download", _capture)
    download.fetch_gencode_gtf("GRCh38", 45, str(tmp_path))
    assert "release_45" in seen["url"]
    assert "gencode.v45" in seen["url"]


def test_fetch_gencode_gtf_url_for_mouse(monkeypatch, tmp_path):
    seen = {}

    def _capture(url, dst):
        seen["url"] = url
        with open(dst, "wb") as f:
            f.write(b"x")
        return dst

    monkeypatch.setattr(download, "_download", _capture)
    download.fetch_gencode_gtf("GRCm39", 35, str(tmp_path))
    assert "release_M35" in seen["url"]
    assert "gencode.vM35" in seen["url"]
