"""Every CATEGORY and BUCKET must have a complete glossary entry, since
README, HTML report, and `tecap explain` all render from these dicts."""

import re

from tecap.constants import (
    BUCKET_DEFINITIONS,
    BUCKETS,
    CATEGORIES,
    GLOSSARY,
    HOW_TO_READ_COMPARE,
    MECHANISM_DEFINITIONS,
    PLOT_CAPTIONS,
    REPORT_INTRO,
)

REQUIRED = ("short", "what", "why")

# A bucket says where the site is, the label says what kind of site it is:
# no adenine band and no bucket names a priming mechanism. Written as a
# regex so the repo-wide grep for the retired phrases stays empty.
RETIRED = re.compile(r"\bmoderate.a\b|\bclassical.a\b|\ba.tract|mispriming", re.IGNORECASE)


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


def test_no_text_names_a_priming_mechanism_from_a_band_or_bucket():
    texts = [v for d in (*MECHANISM_DEFINITIONS.values(), *BUCKET_DEFINITIONS.values())
             for v in d.values()]
    texts += [*PLOT_CAPTIONS.values(), *REPORT_INTRO.values(), HOW_TO_READ_COMPARE]
    texts += [t for row in GLOSSARY for t in row]
    for text in texts:
        assert not RETIRED.search(text), text
