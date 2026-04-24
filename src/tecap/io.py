"""JSON read/write with schema versioning."""

import json

from tecap.constants import SCHEMA_VERSION


def write_json(path, payload):
    """Write a summary dict with schema_version stamped in."""
    out = dict(payload)
    out.setdefault("schema_version", SCHEMA_VERSION)
    with open(path, "w") as f:
        json.dump(out, f, indent=2)


def read_json(path):
    with open(path) as f:
        return json.load(f)
