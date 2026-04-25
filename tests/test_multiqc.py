"""Tests for the MultiQC custom-content adapter."""

import json

from tecap.constants import (
    CAPTURED,
    MECH_A_CORRECT,
    MECH_B_APA,
    MECH_B_ASPECI,
)
from tecap.multiqc import build_mqc_payload, write_mqc_json


def _fake_results():
    return {
        "total": 200,
        "counts": {
            CAPTURED:        80,
            MECH_A_CORRECT:  40,
            MECH_B_APA:      30,
            MECH_B_ASPECI:   20,
        },
    }


def _fake_summary():
    return {
        "pas_fraction":      {MECH_A_CORRECT: 0.75, MECH_B_APA: 0.25},
        "orient_match_frac": 0.98,
    }


def test_build_mqc_payload_shape():
    p = build_mqc_payload("S1", _fake_results(), _fake_summary())
    assert p["plot_type"] == "table"
    assert "S1" in p["data"]
    row = p["data"]["S1"]
    assert row["captured_pct"]      == 40.0
    assert row["mecha_correct_pct"] == 20.0
    assert row["mechb_aspeci_pct"]  == 10.0
    assert row["mecha_pas_pos_frac"]     == 0.75
    assert row["mechb_apa_pas_pos_frac"] == 0.25
    assert row["orient_match_frac"]      == 0.98
    assert row["total_classified"]       == 200


def test_write_mqc_json_roundtrip(tmp_path):
    p = build_mqc_payload("S1", _fake_results(), _fake_summary())
    out = tmp_path / "S1_tecap_mqc.json"
    write_mqc_json(str(out), p)
    loaded = json.loads(out.read_text())
    assert loaded["id"] == "tecap_classify"
    assert loaded["data"]["S1"]["captured_pct"] == 40.0


def test_build_mqc_payload_empty_total():
    res = {"total": 0, "counts": {}}
    p = build_mqc_payload("S0", res, {"pas_fraction": {}, "orient_match_frac": None})
    row = p["data"]["S0"]
    assert row["captured_pct"]       is None
    assert row["mecha_correct_pct"]  is None
    assert row["total_classified"]   == 0
