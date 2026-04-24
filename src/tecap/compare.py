"""Cross-sample comparison dispatcher for classify and basecomp JSON outputs."""

import logging

from tecap.io import read_json
from tecap.plotting import plot_comparison_basecomp, plot_comparison_classify

log = logging.getLogger(__name__)


def run_compare(mode, input_paths, out_dir):
    if mode == "classify":
        all_results = {}
        for p in input_paths:
            d = read_json(p)
            all_results[d["sample"]] = d
        plot_comparison_classify(all_results, out_dir)
    elif mode == "basecomp":
        all_summaries = {}
        for p in input_paths:
            d = read_json(p)
            all_summaries[d["sample"]] = d
        window = next(iter(all_summaries.values())).get("window", 20)
        plot_comparison_basecomp(all_summaries, out_dir, window)
    else:
        raise ValueError(f"Unknown compare mode: {mode!r}")
