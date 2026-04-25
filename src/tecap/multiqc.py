"""MultiQC custom-content adapter. Writes `{sample}_tecap_mqc.json` next to
the main classify JSON. MultiQC auto-detects the `_mqc.json` suffix and
renders a per-sample table."""

import json

from tecap.constants import (
    CAPTURED,
    MECH_A_CORRECT,
    MECH_B_APA,
    MECH_B_ASPECI,
)


def build_mqc_payload(sample, results, summary):
    """Compose a MultiQC custom-content table payload from a classify result.

    `results` is the analyse_bam return dict; `summary` is the JSON summary
    actually written for the user (contains pas_fraction, orient_match_frac).
    """
    total = results["total"]
    counts = results["counts"]

    def pct(cat):
        return (counts.get(cat, 0) / total * 100) if total > 0 else None

    pas_frac = summary.get("pas_fraction") or {}
    orient_match_frac = summary.get("orient_match_frac")

    row = {
        "captured_pct":               pct(CAPTURED),
        "mecha_correct_pct":          pct(MECH_A_CORRECT),
        "mechb_aspeci_pct":           pct(MECH_B_ASPECI),
        "orient_match_frac":          orient_match_frac,
        "mecha_pas_pos_frac":         pas_frac.get(MECH_A_CORRECT),
        "mechb_apa_pas_pos_frac":     pas_frac.get(MECH_B_APA),
        "total_classified":           total,
    }

    return {
        "id": "tecap_classify",
        "section_name": "tecap",
        "description": (
            "Terminal-exon capture diagnostics: bucket fractions, PAS hexamer "
            "split, and BAM/gene-strand orientation sanity check."
        ),
        "plot_type": "table",
        "pconfig": {
            "id": "tecap_classify_table",
            "title": "tecap: terminal-exon capture",
        },
        "headers": {
            "captured_pct": {
                "title": "Captured %",
                "description": "Reads covering >=50% of the terminal exon",
                "min": 0, "max": 100, "suffix": "%",
                "format": "{:,.1f}",
            },
            "mecha_correct_pct": {
                "title": "MechA-correct %",
                "description": "3' end at a polyA site in the TE UTR, but read too short",
                "min": 0, "max": 100, "suffix": "%",
                "format": "{:,.1f}",
            },
            "mechb_aspeci_pct": {
                "title": "MechB-aspecific %",
                "description": "Internal priming upstream of TE, not at any known polyA site",
                "min": 0, "max": 100, "suffix": "%",
                "format": "{:,.1f}",
            },
            "orient_match_frac": {
                "title": "Orient match",
                "description": "Fraction of reads whose BAM strand matches gene strand",
                "min": 0, "max": 1,
                "format": "{:,.3f}",
            },
            "mecha_pas_pos_frac": {
                "title": "MechA PAS+ frac",
                "description": "Fraction of MechA-correct reads at a PAS+ cluster",
                "min": 0, "max": 1,
                "format": "{:,.3f}",
            },
            "mechb_apa_pas_pos_frac": {
                "title": "MechB-APA PAS+ frac",
                "description": "Fraction of MechB-APA reads at a PAS+ cluster",
                "min": 0, "max": 1,
                "format": "{:,.3f}",
            },
            "total_classified": {
                "title": "Reads",
                "description": "Total reads landing in any bucket",
                "format": "{:,d}",
            },
        },
        "data": {sample: row},
    }


def write_mqc_json(path, payload):
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)
