"""Helpers for fetching reference resources (PolyASite atlas, GENCODE GTF).

URLs are stdlib-fetched (no extra deps). SHA256 is computed and logged on
download; users can pin a checksum via --expected-sha256 to verify
reproducibility across runs."""

import hashlib
import logging
import os
import urllib.request

log = logging.getLogger(__name__)


POLYA_URLS = {
    "GRCh38": (
        "https://polyasite.unibas.ch/download/atlas/3.0/"
        "GRCh38.GENCODE_42/atlas.clusters.3.0.GRCh38.GENCODE_42.bed.gz"
    ),
    "GRCm38": (
        "https://polyasite.unibas.ch/download/atlas/2.0/"
        "GRCm38.96/atlas.clusters.2.0.GRCm38.96.bed.gz"
    ),
}


GENCODE_HUMAN_TEMPLATE = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_{ver}/gencode.v{ver}.annotation.gtf.gz"
)
GENCODE_MOUSE_TEMPLATE = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/"
    "release_M{ver}/gencode.vM{ver}.annotation.gtf.gz"
)


def _sha256_of(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            buf = f.read(chunk)
            if not buf:
                break
            h.update(buf)
    return h.hexdigest()


def _download(url, dst):
    log.info("Downloading %s", url)
    tmp = dst + ".part"
    with urllib.request.urlopen(url) as resp, open(tmp, "wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    os.replace(tmp, dst)
    return dst


def fetch_polya_atlas(genome, out_dir, expected_sha256=None, url=None):
    """Fetch the PolyASite 3.0 cluster BED for `genome` into `out_dir`.

    Returns (path, sha256). If `expected_sha256` is given and does not match,
    raises ValueError after writing the file (caller can decide to delete).
    """
    if url is None:
        if genome not in POLYA_URLS:
            raise ValueError(
                f"No built-in URL for genome={genome!r}. "
                f"Pass --url. Known: {sorted(POLYA_URLS)}"
            )
        url = POLYA_URLS[genome]

    os.makedirs(out_dir, exist_ok=True)
    fname = os.path.basename(url)
    dst = os.path.join(out_dir, fname)

    if os.path.exists(dst):
        log.info("Already exists: %s (skipping download)", dst)
    else:
        _download(url, dst)

    sha = _sha256_of(dst)
    log.info("SHA256: %s  %s", sha, dst)

    if expected_sha256 and sha != expected_sha256:
        raise ValueError(
            f"SHA256 mismatch for {dst}: expected {expected_sha256}, got {sha}"
        )
    return dst, sha


def fetch_gencode_gtf(genome, version, out_dir, expected_sha256=None):
    """Fetch a GENCODE GTF for human (GRCh38) or mouse (GRCm39).

    `version` is the integer release number, e.g. 45 for human, 35 for mouse.
    GRCm38 is unsupported here — current GENCODE mouse releases (M26+) are on
    GRCm39. If you need a GRCm38 GTF, fetch GENCODE M25 manually.
    """
    if genome == "GRCh38":
        url = GENCODE_HUMAN_TEMPLATE.format(ver=version)
    elif genome == "GRCm39":
        url = GENCODE_MOUSE_TEMPLATE.format(ver=version)
    else:
        raise ValueError(
            f"GENCODE not configured for genome={genome!r}. "
            f"Supported: GRCh38, GRCm39. (GRCm38 GENCODE is M25-and-earlier; fetch manually.)"
        )

    os.makedirs(out_dir, exist_ok=True)
    fname = os.path.basename(url)
    dst = os.path.join(out_dir, fname)

    if os.path.exists(dst):
        log.info("Already exists: %s (skipping download)", dst)
    else:
        _download(url, dst)

    sha = _sha256_of(dst)
    log.info("SHA256: %s  %s", sha, dst)

    if expected_sha256 and sha != expected_sha256:
        raise ValueError(
            f"SHA256 mismatch for {dst}: expected {expected_sha256}, got {sha}"
        )
    return dst, sha
