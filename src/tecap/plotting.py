"""All matplotlib plotting code for classify + basecomp outputs."""

import logging
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tecap.constants import (
    CATEGORIES, CAPTURED, MECH_A_CORRECT,
    LABELS_SHORT, COLORS, UTR_BIN_LABELS,
    BUCKETS, BUCKET_COLORS,
)

log = logging.getLogger(__name__)


def plot_single(results, sample_name, out_dir, cov_threshold):
    fig, axes = plt.subplots(1, 3, figsize=(20, 7))
    fig.suptitle(f"3' Terminal Exon Capture Analysis -- {sample_name}",
                 fontsize=13, fontweight="bold")

    ax     = axes[0]
    fracs  = [results["fractions"].get(c, 0) * 100 for c in CATEGORIES]
    labels = [LABELS_SHORT[c] for c in CATEGORIES]
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
                alpha=0.6, label="MechA-correct", color=COLORS[MECH_A_CORRECT], density=True)
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
        sc_path = os.path.join(out_dir, f"{sample_name}_mecha_scatter.png")
        plt.savefig(sc_path, dpi=150, bbox_inches="tight")
        plt.close()
        log.info("Saved plot: %s", sc_path)


def plot_comparison_classify(all_results, out_dir):
    samples = list(all_results.keys())
    n       = len(samples)

    fig, axes = plt.subplots(1, 3, figsize=(22, 7))
    fig.suptitle("3' Terminal Exon Capture -- Multi-sample Comparison",
                 fontsize=13, fontweight="bold")

    ax    = axes[0]
    x     = np.arange(len(CATEGORIES))
    width = 0.8 / n
    for i, sample in enumerate(samples):
        fracs = [all_results[sample]["fractions"].get(c, 0) * 100 for c in CATEGORIES]
        ax.bar(x + i * width, fracs, width, label=sample, alpha=0.85)
    ax.set_xticks(x + width * (n - 1) / 2)
    ax.set_xticklabels([LABELS_SHORT[c] for c in CATEGORIES], fontsize=7)
    ax.set_ylabel("% of multi-exon reads")
    ax.set_title("All mechanisms by sample")
    ax.legend(fontsize=9)

    ax     = axes[1]
    x2     = np.arange(len(UTR_BIN_LABELS))
    width2 = 0.8 / n
    sample_colors = plt.cm.Set2(np.linspace(0, 0.7, n))
    for i, sample in enumerate(samples):
        rates = [v * 100 for v in all_results[sample]["utr_bin_mecha_rate"]]
        ax.bar(x2 + i * width2, rates, width2, label=sample,
               color=sample_colors[i], alpha=0.85)
    ax.set_xticks(x2 + width2 * (n - 1) / 2)
    ax.set_xticklabels(UTR_BIN_LABELS, fontsize=9)
    ax.set_xlabel("3' UTR length bin (bp)")
    ax.set_ylabel("MechA-correct rate\n(% of reads in bin)")
    ax.set_title("MechA-correct rate by UTR length bin")
    ax.legend(fontsize=9)

    ax = axes[2]
    for i, sample in enumerate(samples):
        rates = [v * 100 for v in all_results[sample]["utr_bin_capt_rate"]]
        ax.bar(x2 + i * width2, rates, width2, label=sample,
               color=sample_colors[i], alpha=0.85)
    ax.set_xticks(x2 + width2 * (n - 1) / 2)
    ax.set_xticklabels(UTR_BIN_LABELS, fontsize=9)
    ax.set_xlabel("3' UTR length bin (bp)")
    ax.set_ylabel("Captured rate\n(% of reads in bin)")
    ax.set_title("Captured rate by UTR length bin")
    ax.legend(fontsize=9)

    plt.tight_layout()
    out_path = os.path.join(out_dir, "comparison_terminal_exon.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info("Saved plot: %s", out_path)


def plot_histograms_basecomp(summary, sample, out_dir, window):
    bin_edges = np.asarray(summary["bin_edges"])
    bin_mid   = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    fig, axes = plt.subplots(2, 4, figsize=(18, 8), sharex=True, sharey=True)
    fig.suptitle(f"Base composition {window} nt downstream of cleavage "
                 f"(reference, mRNA orientation) -- {sample}",
                 fontsize=13, fontweight="bold")
    axes = axes.flatten()

    for ax, bucket in zip(axes, BUCKETS):
        data = summary["buckets"][bucket]
        n    = data["n"]
        if n == 0:
            ax.set_title(f"{bucket}  (n=0)", fontsize=10)
            ax.axis("off")
            continue
        hist = np.asarray(data["hist"])
        hist_n = hist / hist.sum() if hist.sum() > 0 else hist
        ax.bar(bin_mid, hist_n, width=0.045,
               color=BUCKET_COLORS[bucket], edgecolor="white", linewidth=0.3)
        ax.axvspan(0.30, 0.50, color="grey", alpha=0.12,
                   label="30-50%" if bucket == BUCKETS[0] else None)
        ax.axvline(0.60, color="black", linestyle="--", linewidth=0.8,
                   label=">=60%" if bucket == BUCKETS[0] else None)
        ax.set_title(
            f"{bucket}  (n={n:,})\n"
            f"med={data['median']*100:.1f}%  "
            f">=60%: {data['frac_ge_60']*100:.1f}%  "
            f"30-50%: {data['frac_30_50']*100:.1f}%",
            fontsize=9,
        )
        ax.set_xlim(0, 1)
    for ax in axes[-4:]:
        ax.set_xlabel(f"%A in {window} nt window")
    for ax in (axes[0], axes[4]):
        ax.set_ylabel("Fraction of reads")

    plt.tight_layout()
    out_path = os.path.join(out_dir, f"{sample}_basecomp.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info("Saved plot: %s", out_path)


def plot_comparison_basecomp(all_summaries, out_dir, window):
    samples = list(all_summaries.keys())
    n_samp  = len(samples)
    bin_edges = np.asarray(next(iter(all_summaries.values()))["bin_edges"])
    bin_mid   = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    width     = (bin_edges[1] - bin_edges[0]) * 0.45

    fig, axes = plt.subplots(2, 4, figsize=(20, 9), sharex=True, sharey=True)
    fig.suptitle(f"Base composition {window} nt downstream of cleavage -- cross-sample comparison",
                 fontsize=13, fontweight="bold")
    axes = axes.flatten()

    palette = plt.cm.Set2(np.linspace(0, 0.7, max(n_samp, 1)))
    sample_colors = {s: palette[i] for i, s in enumerate(samples)}

    for ax, bucket in zip(axes, BUCKETS):
        for i, s in enumerate(samples):
            data = all_summaries[s]["buckets"][bucket]
            if data["n"] == 0:
                continue
            hist = np.asarray(data["hist"], dtype=float)
            hist_n = hist / hist.sum() if hist.sum() > 0 else hist
            off = (i - (n_samp - 1) / 2.0) * width
            ax.bar(bin_mid + off, hist_n, width=width,
                   color=sample_colors[s], alpha=0.8,
                   label=f"{s} (n={data['n']:,}, med={data['median']*100:.1f}%)")
        ax.axvspan(0.30, 0.50, color="grey", alpha=0.12)
        ax.axvline(0.60, color="black", linestyle="--", linewidth=0.8)
        ax.set_title(bucket, fontsize=10)
        ax.set_xlim(0, 1)
        ax.legend(fontsize=7, loc="upper right")
    for ax in axes[-4:]:
        ax.set_xlabel("%A in window")
    for ax in (axes[0], axes[4]):
        ax.set_ylabel("Fraction of reads")

    plt.tight_layout()
    out_path = os.path.join(out_dir, "comparison_basecomp.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info("Saved plot: %s", out_path)
