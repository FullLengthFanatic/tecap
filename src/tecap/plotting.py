"""All matplotlib plotting code for classify + basecomp outputs."""

import logging
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from tecap.constants import (
    BUCKET_COLORS,
    BUCKET_DEFINITIONS,
    BUCKETS,
    CAPTURED,
    CATEGORIES,
    COLORS,
    MECH_A_CORRECT,
    MECHANISM_DEFINITIONS,
    PLOT_CAPTIONS,
    UTR_BIN_LABELS,
)

log = logging.getLogger(__name__)

CAPTION_KW = dict(ha="center", va="top", fontsize=8.5, style="italic", color="#444",
                  wrap=True)

# Cross-sample palette: tab10 reordered to skip the bright olive (#bcbd22)
# and the gray (#7f7f7f), neither of which read well on white.
SAMPLE_PALETTE = [
    "#1f77b4",  # blue
    "#ff7f0e",  # orange
    "#2ca02c",  # green
    "#d62728",  # red
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
    "#17becf",  # teal
]

# Distinct from any green so the read-length histogram does not collide with
# the Captured / MechA-correct bar colors.
HIST_MECHA_COLOR = "#d35400"


def _sample_color(i):
    return SAMPLE_PALETTE[i % len(SAMPLE_PALETTE)]


def _add_caption(fig, text):
    fig.text(0.5, -0.01, text, **CAPTION_KW)


def _basecomp_legend_handles():
    return [
        Patch(facecolor="grey", alpha=0.25,
              label="30-50% A: moderate-A priming"),
        Line2D([0], [0], color="black", linestyle="--", linewidth=1.0,
               label=">=60% A: classical A-tract"),
    ]


def plot_single(results, sample_name, out_dir, cov_threshold):
    fig, axes = plt.subplots(1, 3, figsize=(20, 7))
    fig.suptitle(f"3' Terminal Exon Capture Analysis -- {sample_name}",
                 fontsize=13, fontweight="bold")

    ax     = axes[0]
    fracs  = [results["fractions"].get(c, 0) * 100 for c in CATEGORIES]
    labels = [MECHANISM_DEFINITIONS[c]["short"] for c in CATEGORIES]
    colors = [COLORS[c] for c in CATEGORIES]
    bars   = ax.barh(labels, fracs, color=colors, edgecolor="white", linewidth=0.5)
    for bar, frac in zip(bars, fracs):
        if frac > 0.3:
            ax.text(frac + 0.3, bar.get_y() + bar.get_height() / 2,
                    f"{frac:.1f}%", va="center", fontsize=8)
    ax.set_xlabel("% of multi-exon reads")
    ax.set_xlim(0, max(fracs) * 1.20 if max(fracs) > 0 else 1)
    ax.set_title("Mechanism breakdown")
    ax.invert_yaxis()

    ax = axes[1]
    max_rl = 5000
    bins   = np.linspace(0, max_rl, 60)
    if results["capt_read_lengths"]:
        ax.hist([min(v, max_rl) for v in results["capt_read_lengths"]], bins=bins,
                alpha=0.6, label="Captured", color=COLORS[CAPTURED], density=True)
    if results["mecha_read_lengths"]:
        ax.hist([min(v, max_rl) for v in results["mecha_read_lengths"]], bins=bins,
                alpha=0.6, label="MechA-correct", color=HIST_MECHA_COLOR, density=True)
    ax.set_xlabel("Read length (bp, capped 5 kb)")
    ax.set_ylabel("Density")
    ax.set_title("Read length: Captured vs MechA-correct")
    ax.legend(fontsize=8)

    ax      = axes[2]
    x       = np.arange(len(UTR_BIN_LABELS))
    width   = 0.35
    mecha_r = [v * 100 for v in results["utr_bin_mecha_rate"]]
    capt_r  = [v * 100 for v in results["utr_bin_capt_rate"]]
    ax.bar(x - width/2, capt_r,  width, label="Captured",     color=COLORS[CAPTURED],       alpha=0.85)
    ax.bar(x + width/2, mecha_r, width, label="MechA-correct", color=COLORS[MECH_A_CORRECT], alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(UTR_BIN_LABELS, fontsize=9)
    ax.set_xlabel("3' UTR length bin (bp)")
    ax.set_ylabel("% of reads in bin")
    ax.set_title("Captured vs MechA-correct rate\nby UTR length bin")
    ax.legend(fontsize=8)

    plt.tight_layout()
    _add_caption(fig, PLOT_CAPTIONS["single"].format(sample=sample_name))
    out_path = os.path.join(out_dir, f"{sample_name}_terminal_exon.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info("Saved plot: %s", out_path)

    if results["mecha_cov_scatter"]:
        fig2, ax2 = plt.subplots(figsize=(7, 5))
        rl  = [v[0] for v in results["mecha_cov_scatter"]]
        cov = [v[1] for v in results["mecha_cov_scatter"]]
        hb = ax2.hexbin(rl, cov, gridsize=50, cmap="YlOrRd",
                        extent=(0, 5000, 0, cov_threshold), mincnt=1)
        plt.colorbar(hb, ax=ax2, label="Read count")
        ax2.axhline(cov_threshold, color="black", linestyle="--",
                    linewidth=1, label=f"Threshold ({cov_threshold})")
        ax2.set_xlabel("Read length (bp)")
        ax2.set_ylabel("Coverage fraction of terminal exon")
        ax2.set_title(f"MechA-correct: read length vs TE coverage\n{sample_name}")
        ax2.legend(fontsize=8)
        plt.tight_layout()
        _add_caption(fig2, PLOT_CAPTIONS["mecha_scatter"])
        sc_path = os.path.join(out_dir, f"{sample_name}_mecha_scatter.png")
        plt.savefig(sc_path, dpi=150, bbox_inches="tight")
        plt.close()
        log.info("Saved plot: %s", sc_path)


def plot_comparison_classify(all_results, out_dir):
    samples = list(all_results.keys())
    n       = len(samples)

    fig_w = 22 + 2 * max(0, n - 2)
    fig, axes = plt.subplots(1, 3, figsize=(fig_w, 8))
    fig.suptitle("3' Terminal Exon Capture -- Multi-sample Comparison",
                 fontsize=13, fontweight="bold")

    # Panel 0: horizontal grouped bars, one row per category, samples within.
    ax            = axes[0]
    y             = np.arange(len(CATEGORIES))
    bar_h         = 0.8 / n
    long_labels   = [MECHANISM_DEFINITIONS[c]["short"] for c in CATEGORIES]
    for i, sample in enumerate(samples):
        fracs  = [all_results[sample]["fractions"].get(c, 0) * 100 for c in CATEGORIES]
        offset = (i - (n - 1) / 2.0) * bar_h
        bars   = ax.barh(y + offset, fracs, height=bar_h,
                         label=sample, color=_sample_color(i),
                         edgecolor="white", linewidth=0.4)
        for bar, frac in zip(bars, fracs):
            if frac > 0.5:
                ax.text(frac + 0.3, bar.get_y() + bar.get_height() / 2,
                        f"{frac:.1f}%", va="center", fontsize=7)
    ax.set_yticks(y)
    ax.set_yticklabels(long_labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("% of multi-exon reads")
    ax.set_title("All mechanisms by sample")
    ax.legend(fontsize=9, loc="lower right")

    # Panels 1 + 2: line+marker per sample over UTR length bins. Lines stay
    # readable at any N and don't hide samples whose rate is zero in a bin.
    ax = axes[1]
    x2 = np.arange(len(UTR_BIN_LABELS))
    for i, sample in enumerate(samples):
        rates = [v * 100 for v in all_results[sample]["utr_bin_mecha_rate"]]
        ax.plot(x2, rates, marker="o", linewidth=2, markersize=6,
                label=sample, color=_sample_color(i), alpha=0.9)
    ax.set_xticks(x2)
    ax.set_xticklabels(UTR_BIN_LABELS, fontsize=9)
    ax.set_xlabel("3' UTR length bin (bp)")
    ax.set_ylabel("MechA-correct rate (% of reads in bin)")
    ax.set_title("MechA-correct rate by UTR length bin")
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    ax.legend(fontsize=9)

    ax = axes[2]
    for i, sample in enumerate(samples):
        rates = [v * 100 for v in all_results[sample]["utr_bin_capt_rate"]]
        ax.plot(x2, rates, marker="o", linewidth=2, markersize=6,
                label=sample, color=_sample_color(i), alpha=0.9)
    ax.set_xticks(x2)
    ax.set_xticklabels(UTR_BIN_LABELS, fontsize=9)
    ax.set_xlabel("3' UTR length bin (bp)")
    ax.set_ylabel("Captured rate (% of reads in bin)")
    ax.set_title("Captured rate by UTR length bin")
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    ax.legend(fontsize=9)

    plt.tight_layout()
    _add_caption(fig, PLOT_CAPTIONS["comparison_classify"])
    out_path = os.path.join(out_dir, "comparison_terminal_exon.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info("Saved plot: %s", out_path)


def plot_histograms_basecomp(summary, sample, out_dir, window):
    bin_edges = np.asarray(summary["bin_edges"])
    bin_mid   = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    fig, axes = plt.subplots(2, 4, figsize=(22, 11), sharex=True, sharey=True)
    fig.suptitle(f"Base composition {window} nt downstream of cleavage "
                 f"(reference, mRNA orientation) -- {sample}",
                 fontsize=15, fontweight="bold")
    axes = axes.flatten()

    for ax, bucket in zip(axes, BUCKETS):
        data = summary["buckets"][bucket]
        n    = data["n"]
        if n == 0:
            ax.set_title(f"{bucket}  (n=0)", fontsize=11)
            ax.axis("off")
            continue
        hist = np.asarray(data["hist"])
        hist_n = hist / hist.sum() if hist.sum() > 0 else hist
        ax.bar(bin_mid, hist_n, width=0.045,
               color=BUCKET_COLORS[bucket], edgecolor="white", linewidth=0.3)
        ax.axvspan(0.30, 0.50, color="grey", alpha=0.18)
        ax.axvline(0.60, color="black", linestyle="--", linewidth=0.8)
        ax.set_title(
            f"{bucket}  (n={n:,})\n"
            f"med={data['median']*100:.1f}%  "
            f">=60%: {data['frac_ge_60']*100:.1f}%  "
            f"30-50%: {data['frac_30_50']*100:.1f}%",
            fontsize=11, fontweight="bold",
        )
        ax.tick_params(axis="both", labelsize=10)
        ax.set_xlim(0, 1)
    for ax in axes[-4:]:
        ax.set_xlabel(f"%A in {window} nt window", fontsize=11)
    for ax in (axes[0], axes[4]):
        ax.set_ylabel("Fraction of reads", fontsize=11)

    glossary = "  |  ".join(
        f"{b}: {BUCKET_DEFINITIONS[b]['why']}" for b in BUCKETS
    )
    fig.text(0.5, -0.04, glossary, ha="center", va="top",
             fontsize=8.5, color="#444", wrap=True)

    fig.legend(handles=_basecomp_legend_handles(),
               loc="upper right", bbox_to_anchor=(0.995, 0.985), fontsize=10,
               framealpha=0.9)

    plt.tight_layout(rect=(0, 0, 1, 0.96))
    _add_caption(fig, PLOT_CAPTIONS["basecomp"].format(window=window))
    out_path = os.path.join(out_dir, f"{sample}_basecomp.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info("Saved plot: %s", out_path)


def plot_comparison_basecomp(all_summaries, out_dir, window):
    samples = list(all_summaries.keys())
    bin_edges = np.asarray(next(iter(all_summaries.values()))["bin_edges"])
    bin_mid   = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    fig, axes = plt.subplots(2, 4, figsize=(24, 12), sharex=True, sharey=True)
    fig.suptitle(f"Base composition {window} nt downstream of cleavage -- cross-sample comparison",
                 fontsize=15, fontweight="bold")
    axes = axes.flatten()

    sample_colors = {s: _sample_color(i) for i, s in enumerate(samples)}

    legend_handles = {}
    for ax, bucket in zip(axes, BUCKETS):
        for s in samples:
            data = all_summaries[s]["buckets"][bucket]
            if data["n"] == 0:
                continue
            hist = np.asarray(data["hist"], dtype=float)
            hist_n = hist / hist.sum() if hist.sum() > 0 else hist
            label = f"{s} (n={data['n']:,}, med={data['median']*100:.1f}%)"
            line, = ax.step(bin_mid, hist_n, where="mid", linewidth=1.8,
                            color=sample_colors[s], alpha=0.9,
                            label=label)
            ax.fill_between(bin_mid, hist_n, step="mid",
                            color=sample_colors[s], alpha=0.10)
            legend_handles.setdefault(s, line)
        ax.axvspan(0.30, 0.50, color="grey", alpha=0.18)
        ax.axvline(0.60, color="black", linestyle="--", linewidth=0.8)
        ax.set_title(bucket, fontsize=12, fontweight="bold")
        ax.set_xlim(0, 1)
        ax.tick_params(axis="both", labelsize=10)
    for ax in axes[-4:]:
        ax.set_xlabel(f"%A in {window} nt window", fontsize=11)
    for ax in (axes[0], axes[4]):
        ax.set_ylabel("Fraction of reads", fontsize=11)

    if legend_handles:
        fig.legend(handles=list(legend_handles.values()),
                   labels=list(legend_handles.keys()),
                   loc="upper center", bbox_to_anchor=(0.5, 0.965),
                   ncol=min(len(legend_handles), 4), fontsize=10,
                   framealpha=0.9)

    # Bucket glossary at the bottom of the figure so subplot titles stay
    # uncluttered but readers still see what each bucket means.
    glossary = "  |  ".join(
        f"{b}: {BUCKET_DEFINITIONS[b]['why']}" for b in BUCKETS
    )
    fig.text(0.5, -0.04, glossary, ha="center", va="top",
             fontsize=8.5, color="#444", wrap=True)

    fig.legend(handles=_basecomp_legend_handles(),
               loc="upper right", bbox_to_anchor=(0.995, 0.965), fontsize=9,
               framealpha=0.9)

    plt.tight_layout(rect=(0, 0, 1, 0.92))
    _add_caption(fig, PLOT_CAPTIONS["comparison_basecomp"].format(window=window))
    out_path = os.path.join(out_dir, "comparison_basecomp.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info("Saved plot: %s", out_path)


