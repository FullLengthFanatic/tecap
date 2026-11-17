"""HTML report tests against synthetic classify+basecomp JSON fixtures."""

from html.parser import HTMLParser

import pytest

from tecap.constants import BUCKETS, CATEGORIES, MECHANISM_DEFINITIONS, SCHEMA_VERSION
from tecap.io import write_json
from tecap.report import build_compare_report, build_single_report


def _classify_summary(sample, total=10000):
    return {
        "schema_version": SCHEMA_VERSION,
        "sample": sample,
        "platform": "cdna-pacbio",
        "total": total,
        "counts": {c: total // len(CATEGORIES) for c in CATEGORIES},
        "fractions": {c: 1.0 / len(CATEGORIES) for c in CATEGORIES},
        "polya_slop": 25,
        "polya_types": "TE,DS,DI,EX",
        "pas_counts": {
            CATEGORIES[1]: {"pas": 700, "nopas": 300},
            CATEGORIES[5]: {"pas": 200, "nopas": 100},
        },
        "pas_fraction": {CATEGORIES[1]: 0.7, CATEGORIES[5]: 0.667},
        "orient_match_frac": 0.972,
        "utr_bin_labels":         ["0-500", "500-1k", "1k-2k", "2k-5k", ">5k"],
        "utr_bin_total":          [100, 200, 300, 400, 500],
        "utr_bin_mecha_correct":  [10, 20, 30, 40, 50],
        "utr_bin_captured":       [50, 100, 150, 200, 250],
        "utr_bin_mecha_rate":     [0.10, 0.10, 0.10, 0.10, 0.10],
        "utr_bin_capt_rate":      [0.50, 0.50, 0.50, 0.50, 0.50],
        "cov_median_by_cat":      {},
        "read_len_median_mecha":  900,
        "read_len_median_capt":   2400,
    }


def _basecomp_summary(sample):
    bin_edges = [i / 20 for i in range(21)]
    return {
        "schema_version": SCHEMA_VERSION,
        "sample": sample,
        "platform": "cdna-pacbio",
        "window": 20,
        "polya_slop": 25,
        "polya_types": "TE,DS,DI,EX",
        "bin_edges": bin_edges,
        "buckets": {
            "Captured":         {"n": 5000, "median": 0.18, "frac_ge_60": 0.05, "frac_30_50": 0.12, "hist": [0]*20},
            "MechA_PAS+":       {"n": 1000, "median": 0.20, "frac_ge_60": 0.06, "frac_30_50": 0.15, "hist": [0]*20},
            "MechA_PAS-":       {"n":  800, "median": 0.40, "frac_ge_60": 0.20, "frac_30_50": 0.50, "hist": [0]*20},
            "MechB_APA_PAS+":   {"n":  300, "median": 0.22, "frac_ge_60": 0.07, "frac_30_50": 0.16, "hist": [0]*20},
            "MechB_APA_PAS-":   {"n":  150, "median": 0.45, "frac_ge_60": 0.30, "frac_30_50": 0.45, "hist": [0]*20},
            "MechB_exon":       {"n":  400, "median": 0.55, "frac_ge_60": 0.40, "frac_30_50": 0.40, "hist": [0]*20},
            "MechB_aspecific":  {"n":  200, "median": 0.62, "frac_ge_60": 0.55, "frac_30_50": 0.30, "hist": [0]*20},
            "IP_TE_CDS":        {"n":  100, "median": 0.50, "frac_ge_60": 0.35, "frac_30_50": 0.40, "hist": [0]*20},
        },
    }


# n_templated, n_nontemplated, n_ambiguous per bucket; n_site is their sum.
_LABEL_COUNTS = {
    "Captured":        (300, 4200, 400),
    "MechA_PAS+":      (60, 820, 90),
    "MechA_PAS-":      (240, 380, 150),
    "MechB_APA_PAS+":  (20, 250, 25),
    "MechB_APA_PAS-":  (55, 60, 30),
    "MechB_exon":      (180, 120, 90),
    "MechB_aspecific": (120, 30, 45),
    "IP_TE_CDS":       (0, 0, 0),
}


def _labelled_basecomp_summary(sample):
    bc = _basecomp_summary(sample)
    for b, (t, nt, amb) in _LABEL_COUNTS.items():
        bc["buckets"][b].update(
            n_site=t + nt + amb, n_templated=t, n_nontemplated=nt, n_ambiguous=amb)
    return bc


def _null_label_basecomp_summary(sample):
    # What basecomp writes when --window is not 20: label keys present, null.
    bc = _basecomp_summary(sample)
    bc["window"] = 30
    for d in bc["buckets"].values():
        d.update(n_site=None, n_templated=None, n_nontemplated=None, n_ambiguous=None)
    return bc


_LABEL_HEADER = ["Bucket", "Reads with a site", "Templated", "Non-templated", "Ambiguous"]

_LABEL_CAPTION = (
    "A bucket says where the site is. The label says what kind of site it is: "
    "templated (8 or more adenines in the 20 nt of reference downstream of the "
    "site), non-templated (4 or fewer) or ambiguous (5 to 7). Counts cover "
    "reads whose site was observed under the site anchor of this run."
)


class _TableCells(HTMLParser):
    """Every <table> as a list of rows, each row a list of cell texts."""

    def __init__(self):
        super().__init__()
        self.tables = []
        self._cell = None

    def handle_starttag(self, tag, attrs):
        if tag == "table":
            self.tables.append([])
        elif tag == "tr":
            self.tables[-1].append([])
        elif tag in ("td", "th"):
            self._cell = []

    def handle_endtag(self, tag):
        if tag in ("td", "th"):
            self.tables[-1][-1].append("".join(self._cell).strip())
            self._cell = None

    def handle_data(self, data):
        if self._cell is not None:
            self._cell.append(data)


def _label_tables(html):
    p = _TableCells()
    p.feed(html)
    return [t for t in p.tables if t and t[0] == _LABEL_HEADER]


def _expected_label_rows():
    rows = []
    for b in BUCKETS:
        t, nt, amb = _LABEL_COUNTS[b]
        rows.append([b, f"{t + nt + amb:,}", f"{t:,}", f"{nt:,}", f"{amb:,}"])
    return rows


class _ImgCounter(HTMLParser):
    def __init__(self):
        super().__init__()
        self.imgs = 0
        self.tables = 0

    def handle_starttag(self, tag, attrs):
        if tag == "img":
            self.imgs += 1
        if tag == "table":
            self.tables += 1


def _count(html):
    p = _ImgCounter()
    p.feed(html)
    return p


def test_single_report_contains_all_mechanisms(tmp_path):
    cj = tmp_path / "S1_terminal_exon.json"
    write_json(str(cj), _classify_summary("S1"))
    html = build_single_report(str(cj))

    for cat in CATEGORIES:
        short = MECHANISM_DEFINITIONS[cat]["short"]
        assert short in html, f"missing short label for {cat}"

    assert "S1" in html
    assert "Executive summary" in html
    assert "Mechanism legend" in html
    assert "schema 0.2" in html.lower() or 'schema 0.2' in html

    counts = _count(html)
    assert counts.tables >= 3


def test_single_report_with_basecomp(tmp_path):
    cj = tmp_path / "S1_terminal_exon.json"
    bj = tmp_path / "S1_basecomp.json"
    write_json(str(cj), _classify_summary("S1"))
    write_json(str(bj), _basecomp_summary("S1"))

    html = build_single_report(str(cj), str(bj))
    assert "Base composition" in html
    assert "Median %A" in html
    assert "30-50% A" in html or "30-50%" in html


def test_single_report_bucket_by_label_table(tmp_path):
    cj = tmp_path / "S1_terminal_exon.json"
    bj = tmp_path / "S1_basecomp.json"
    write_json(str(cj), _classify_summary("S1"))
    write_json(str(bj), _labelled_basecomp_summary("S1"))

    html = build_single_report(str(cj), str(bj))
    tables = _label_tables(html)
    assert len(tables) == 1
    assert tables[0][1:] == _expected_label_rows()
    assert _LABEL_CAPTION in html


@pytest.mark.parametrize("make_bc", [_basecomp_summary, _null_label_basecomp_summary],
                         ids=["no_label_keys", "null_label_keys"])
def test_single_report_without_labels_has_no_label_table(tmp_path, make_bc):
    cj = tmp_path / "S1_terminal_exon.json"
    bj = tmp_path / "S1_basecomp.json"
    write_json(str(cj), _classify_summary("S1"))
    write_json(str(bj), make_bc("S1"))

    html = build_single_report(str(cj), str(bj))
    assert "Median %A" in html
    assert _label_tables(html) == []
    assert _LABEL_CAPTION not in html


def test_compare_report_renders(tmp_path):
    paths = []
    for s in ("S1", "S2"):
        p = tmp_path / f"{s}_terminal_exon.json"
        write_json(str(p), _classify_summary(s))
        paths.append(str(p))

    html = build_compare_report(paths, out_dir=str(tmp_path))
    assert "comparison report" in html.lower()
    assert "S1" in html and "S2" in html
    assert (tmp_path / "comparison_terminal_exon.png").exists()
    counts = _count(html)
    assert counts.imgs >= 1


def test_compare_report_with_basecomp(tmp_path):
    cpaths, bpaths = [], []
    for s in ("S1", "S2"):
        cp = tmp_path / f"{s}_terminal_exon.json"
        bp = tmp_path / f"{s}_basecomp.json"
        write_json(str(cp), _classify_summary(s))
        write_json(str(bp), _basecomp_summary(s))
        cpaths.append(str(cp))
        bpaths.append(str(bp))

    html = build_compare_report(cpaths, bpaths, out_dir=str(tmp_path))
    assert "Base composition" in html
    assert (tmp_path / "comparison_basecomp.png").exists()
    assert _label_tables(html) == []
    assert _LABEL_CAPTION not in html


def test_compare_report_bucket_by_label_table(tmp_path):
    # S1 carries labels, S2 is an older JSON without them: one table, for S1.
    cpaths, bpaths = [], []
    for s, make_bc in (("S1", _labelled_basecomp_summary), ("S2", _basecomp_summary)):
        cp = tmp_path / f"{s}_terminal_exon.json"
        bp = tmp_path / f"{s}_basecomp.json"
        write_json(str(cp), _classify_summary(s))
        write_json(str(bp), make_bc(s))
        cpaths.append(str(cp))
        bpaths.append(str(bp))

    html = build_compare_report(cpaths, bpaths, out_dir=str(tmp_path))
    tables = _label_tables(html)
    assert len(tables) == 1
    assert tables[0][1:] == _expected_label_rows()
    assert html.count(_LABEL_CAPTION) == 1


def test_single_report_has_intro_and_glossary(tmp_path):
    cj = tmp_path / "S1_terminal_exon.json"
    write_json(str(cj), _classify_summary("S1"))
    html = build_single_report(str(cj))

    assert "class='intro'" in html
    assert "summarises" in html and "terminal exon capture diagnostics" in html
    assert "How to read this report" not in html
    assert "<h2>Glossary</h2>" in html
    assert "Untranslated region" in html
    assert "Polyadenylation signal" in html


def test_compare_report_has_intro_howto_and_glossary(tmp_path):
    cpaths = []
    for s in ("S1", "S2"):
        p = tmp_path / f"{s}_terminal_exon.json"
        write_json(str(p), _classify_summary(s))
        cpaths.append(str(p))

    html = build_compare_report(cpaths, out_dir=str(tmp_path))

    assert "class='intro'" in html
    assert "compares" in html and "terminal exon capture diagnostics" in html
    assert "class='howto'" in html
    assert "How to read this report" in html
    assert "Each sample is assigned one colour" in html
    assert "<h2>Glossary</h2>" in html
    assert "Polyadenylation" in html
