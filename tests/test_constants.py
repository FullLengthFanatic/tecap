"""Every CATEGORY and BUCKET must have a complete glossary entry, since
README, HTML report, and `tecap explain` all render from these dicts."""

from tecap.constants import (
    BUCKET_DEFINITIONS,
    BUCKETS,
    CATEGORIES,
    MECHANISM_DEFINITIONS,
    PLOT_CAPTIONS,
)

REQUIRED = ("short", "what", "why")


def test_every_category_has_definition():
    for cat in CATEGORIES:
        assert cat in MECHANISM_DEFINITIONS, f"missing mechanism: {cat!r}"
        d = MECHANISM_DEFINITIONS[cat]
        for f in REQUIRED:
            assert d.get(f), f"{cat}: missing '{f}'"


def test_no_extra_mechanism_keys():
    extra = set(MECHANISM_DEFINITIONS) - set(CATEGORIES)
    assert not extra, f"unexpected mechanisms: {extra}"


def test_every_bucket_has_definition():
    for b in BUCKETS:
        assert b in BUCKET_DEFINITIONS, f"missing bucket: {b!r}"
        d = BUCKET_DEFINITIONS[b]
        for f in REQUIRED:
            assert d.get(f), f"{b}: missing '{f}'"


def test_plot_captions_present():
    expected_keys = {"single", "mecha_scatter", "comparison_classify",
                     "basecomp", "comparison_basecomp"}
    assert expected_keys <= set(PLOT_CAPTIONS)
    assert "{sample}" in PLOT_CAPTIONS["single"]
    assert "{window}" in PLOT_CAPTIONS["basecomp"]
