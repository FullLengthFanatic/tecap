"""HTML report tests against synthetic classify+basecomp JSON fixtures."""

from html.parser import HTMLParser

from tecap.constants import CATEGORIES, MECHANISM_DEFINITIONS, SCHEMA_VERSION
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
    assert "schema 0.1" in html.lower() or 'schema 0.1' in html

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
