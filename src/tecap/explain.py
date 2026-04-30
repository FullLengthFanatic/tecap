"""Pretty-print MECHANISM_DEFINITIONS / BUCKET_DEFINITIONS for the
`tecap explain` subcommand."""

import json

from tecap.constants import (
    BUCKET_DEFINITIONS,
    BUCKETS,
    CATEGORIES,
    MECHANISM_DEFINITIONS,
)


def _match(query, key, definition):
    q = query.lower()
    return q == key.lower() or q == definition["short"].lower()


def _filter(query, mapping, order):
    if query is None:
        return [(k, mapping[k]) for k in order]
    hits = [(k, mapping[k]) for k in order if _match(query, k, mapping[k])]
    return hits


def _format_text(items, header):
    if not items:
        return ""
    lines = [header, "=" * len(header), ""]
    for key, d in items:
        lines.append(f"{d['short']}  ({key})")
        lines.append(f"  what: {d['what']}")
        lines.append(f"  why : {d['why']}")
        lines.append("")
    return "\n".join(lines)


def render(scope="all", mechanism=None, fmt="text"):
    classify_items = (_filter(mechanism, MECHANISM_DEFINITIONS, CATEGORIES)
                      if scope in ("all", "classify") else [])
    basecomp_items = (_filter(mechanism, BUCKET_DEFINITIONS, BUCKETS)
                      if scope in ("all", "basecomp") else [])

    if mechanism is not None and not classify_items and not basecomp_items:
        raise SystemExit(
            f"No definition found for {mechanism!r}. "
            f"Try `tecap explain` with no --mechanism to list all."
        )

    if fmt == "json":
        payload = {}
        if classify_items:
            payload["classify"] = {k: d for k, d in classify_items}
        if basecomp_items:
            payload["basecomp"] = {k: d for k, d in basecomp_items}
        return json.dumps(payload, indent=2)

    chunks = []
    if classify_items:
        chunks.append(_format_text(classify_items, "Classify mechanisms"))
    if basecomp_items:
        chunks.append(_format_text(basecomp_items, "Basecomp buckets"))
    return "\n".join(c for c in chunks if c)
