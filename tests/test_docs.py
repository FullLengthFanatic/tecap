"""README and CHANGELOG sentences on the site anchor, checked against the CLI
and the benchmark definition in tests/data/definition.json."""

import json
from pathlib import Path

from tecap.cli import _build_parser

ROOT = Path(__file__).parent.parent
LONG_READ = json.loads((ROOT / "tests" / "data" / "definition.json").read_text())["site"]["long_read"]


def _block(name, sep, start):
    """The one block of file `name`, split on `sep`, that starts with `start`, whitespace collapsed."""
    (block,) = [b for b in (ROOT / name).read_text().split(sep) if b.startswith(start)]
    return " ".join(block.split())


def _readme_site_anchor():
    return _block("README.md", "\n\n", "Both JSONs record the site anchor.")


def _changelog_site_anchor():
    return _block("CHANGELOG.md", "\n- ", "`--site-anchor {alignment_end,a_run_start}`")


def test_docs_name_the_long_read_filters_tecap_leaves_out():
    for argv in (["classify", "--bam", "x"], ["basecomp", "--bam", "x", "--fasta", "x"]):
        default = _build_parser().parse_args(argv).min_mapq
        for text in (_readme_site_anchor(), _changelog_site_anchor()):
            assert f"does not require MAPQ {LONG_READ['min_mapq']} (`--min-mapq` defaults to {default})" in text
            assert f"{LONG_READ['min_after']['BD67']} nt of bead sequence after the run" in text
            assert "stage 30 applies both before it runs tecap." in text
