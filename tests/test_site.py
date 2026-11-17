"""Read site rule against the shared synthetic vectors (plan Appendix A) and the
benchmark definition in tests/data/definition.json."""

import json
import re
from collections import namedtuple
from pathlib import Path

import pysam
import pytest

from tecap.constants import (
    DEFINITION_SOURCE,
    LABEL_WINDOW,
    NONTEMPLATED_MAX,
    SITE_RUN_MIN,
    SITE_RUN_TOLERANCE,
    TEMPLATED_MIN,
)
from tecap.site import label_of, read_site, terminal_a_run

CONTIG = "chrS"
CONTIG_LEN = 2000
HEADER = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}]})

# Long-read rows of Appendix A, plus strand. Contig chrS is all C except the
# planted bases (start, bases). Sequences use the appendix notation: "C100"
# is 100 C's. site None: no observed site and no label.
Vector = namedtuple("Vector", "planted pos cigar seq site aligned_a clip_a a_window label")
L_VECTORS = {
    "L1": Vector(None, 900, "100M63S", "C100 A30 G33", 1000, 0, 30, 0, "nontemplated"),
    "L2": Vector((1000, "A10"), 900, "110M63S", "C100 A10 A30 G33", 1000, 10, 30, 10, "templated"),
    # at anchor alignment_end: site 1000, templated
    "L3": Vector((1000, "A10"), 900, "100M", "C100", None, None, None, None, None),
    "L4_5": Vector((1002, "A5"), 900, "100M63S", "C100 A30 G33", 1000, 0, 30, 5, "ambiguous"),
    "L4_6": Vector((1002, "A6"), 900, "100M63S", "C100 A30 G33", 1000, 0, 30, 6, "ambiguous"),
    "L4_7": Vector((1002, "A7"), 900, "100M63S", "C100 A30 G33", 1000, 0, 30, 7, "ambiguous"),
    "L5": Vector(None, 900, "100M63H", "C100", None, None, None, None, None),
    "L6": Vector(None, 900, "100M63S", "C100 A5 G A24 G33", 1000, 0, 29, 0, "nontemplated"),
    "L7": Vector((990, "AAAAGAAAAA"), 900, "100M53S", "C90 AAAAGAAAAA A20 G33", 990, 9, 20, 9, "templated"),
    "L8": Vector(None, 900, "100M63S", "C100 A5 G A3 G A21 G32", None, None, None, None, None),
    "L9": Vector((990, "A15"), 900, "105M", "C90 A15", 990, 15, 0, 15, "templated"),
    "L10": Vector(None, 900, "100M42S", "C100 A9 G33", None, None, None, None, None),
    "L11": Vector(None, 900, "100M64S", "C100 G A30 G33", None, None, None, None, None),
    # run found, window of 10 nt: no label, not counted as a read with a site
    "L1_at_contig_end": Vector(None, 1890, "100M10S", "C100 A10", 1990, 0, 10, None, None),
}

_COMP = str.maketrans("ACGT", "TGCA")


def expand(notation):
    out = []
    for tok in notation.split():
        m = re.fullmatch(r"([ACGT])(\d+)", tok)
        out.append(m.group(1) * int(m.group(2)) if m else tok)
    return "".join(out)


def segment(pos, cigar, seq, flag=0):
    a = pysam.AlignedSegment(HEADER)
    a.query_name = "r"
    a.flag = flag
    a.reference_id = 0
    a.reference_start = pos
    a.mapping_quality = 60
    a.cigarstring = cigar
    a.query_sequence = seq
    return a


def build(v, minus):
    """(read, expected site) for a row; minus gives its mirror on the reverse
    strand: reverse complement, reversed CIGAR, flag 16, pos' = 2000 -
    reference_end, site' = 2000 - site."""
    read = segment(v.pos, v.cigar, expand(v.seq))
    if not minus:
        return read, v.site
    cigar = "".join(reversed(re.findall(r"\d+[MIDNSHP=X]", v.cigar)))
    mirror = segment(CONTIG_LEN - read.reference_end, cigar,
                     expand(v.seq).translate(_COMP)[::-1], flag=16)
    return mirror, None if v.site is None else CONTIG_LEN - v.site


def alignment_end(read, minus):
    return read.reference_start if minus else read.reference_end


@pytest.mark.parametrize("name", L_VECTORS)
@pytest.mark.parametrize("minus", [False, True])
def test_a_run_start(name, minus):
    v = L_VECTORS[name]
    read, want_site = build(v, minus)
    site, observed = read_site(read, "-" if minus else "+", "a_run_start")
    assert observed == (v.site is not None)
    assert site == (want_site if observed else alignment_end(read, minus))


@pytest.mark.parametrize("name", [k for k, v in L_VECTORS.items() if v.site is not None])
@pytest.mark.parametrize("minus", [False, True])
def test_run_parts(name, minus):
    v = L_VECTORS[name]
    read, want_site = build(v, minus)
    assert terminal_a_run(read, minus) == (want_site, v.aligned_a, v.clip_a)


def test_run_parts_of_internal_kept():
    assert terminal_a_run(build(L_VECTORS["L2"], False)[0], False) == (1000, 10, 30)


def test_alignment_end_anchor_misses_the_templated_adenines():
    read, _ = build(L_VECTORS["L2"], False)
    assert read_site(read, "+", "alignment_end") == (1010, True)


def test_read_antisense_to_the_gene_falls_back():          # Review Focus 1
    read, _ = build(L_VECTORS["L1"], False)                 # forward read, poly(A) on the right
    assert read_site(read, "-", "a_run_start") == (read.reference_start, False)


def test_run_covering_the_whole_alignment():                # Review Focus 3
    read = segment(pos=990, cigar="15M", seq="A" * 15)
    assert read_site(read, "+", "a_run_start") == (990, True)


def test_labels():
    assert [label_of(k) for k in (0, 4, 5, 7, 8, 20)] == [
        "nontemplated", "nontemplated", "ambiguous", "ambiguous", "templated", "templated"]


def test_constants_equal_the_benchmark_definition():
    d = json.loads((Path(__file__).parent / "data" / "definition.json").read_text())
    assert (LABEL_WINDOW, TEMPLATED_MIN, NONTEMPLATED_MAX) == (
        d["label"]["window"], d["label"]["templated_min"], d["label"]["nontemplated_max"])
    assert (SITE_RUN_MIN, SITE_RUN_TOLERANCE) == (
        d["site"]["long_read"]["run_min"], d["site"]["long_read"]["run_tolerance"])
    assert DEFINITION_SOURCE == f'{d["source"]} {d["version"]}'
