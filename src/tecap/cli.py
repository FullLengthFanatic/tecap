"""Command-line entry point with subcommands: classify, basecomp, compare."""

import argparse
import logging
import os
import sys

import numpy as np

from tecap import __version__
from tecap.constants import CATEGORIES, MECH_A_CORRECT, MECH_B_APA, UTR_BIN_LABELS

PLATFORMS_OK = ("cdna-pacbio", "cdna-ont")
PLATFORMS_REFUSED = {
    "drna-ont": (
        "Direct-RNA (dRNA) sequencing is not supported. dRNA has no reverse "
        "transcription and therefore no oligo-dT priming artifact to diagnose. "
        "Reads are also sequenced 3' to 5', so terminal-exon coverage geometry "
        "does not map onto the classification buckets used here."
    ),
}


def _setup_logging(verbose):
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def _check_platform(platform):
    if platform in PLATFORMS_REFUSED:
        sys.stderr.write(f"ERROR: {PLATFORMS_REFUSED[platform]}\n")
        sys.exit(1)
    if platform not in PLATFORMS_OK:
        sys.stderr.write(
            f"ERROR: unknown --platform {platform!r}. "
            f"Choose one of: {', '.join(PLATFORMS_OK)}\n"
        )
        sys.exit(1)


def _add_common_indexing_args(p):
    p.add_argument("--polya-slop",  type=int, default=25,
                   help="+-bp around each polyA cluster (default 25).")
    p.add_argument("--polya-types", default="TE,DS,DI,EX",
                   help="Comma-separated PolyASite cluster types to keep "
                        "(default TE,DS,DI,EX). Use 'all' to disable filtering.")
    p.add_argument("--min-mapq",           type=int,   default=0)
    p.add_argument("--coverage-threshold", type=float, default=0.5)
    p.add_argument("--platform", default="cdna-pacbio",
                   help="Platform: cdna-pacbio or cdna-ont. drna-ont is refused.")
    p.add_argument("--threads", type=int, default=1,
                   help="Number of worker processes (contig-level parallelism).")
    p.add_argument("--verbose", action="store_true")


def _build_parser():
    parser = argparse.ArgumentParser(
        prog="tecap",
        description="3' terminal exon capture diagnostics for long-read scRNA-seq",
    )
    parser.add_argument("--version", action="version", version=f"tecap {__version__}")
    sub = parser.add_subparsers(dest="cmd", required=True)

    # classify
    pc = sub.add_parser("classify", help="Classify reads into mechanism buckets")
    pc.add_argument("--bam", required=True)
    pc.add_argument("--gtf", required=True)
    pc.add_argument("--polya-sites", required=True)
    pc.add_argument("--sample", default="sample")
    pc.add_argument("--out-dir", default=".")
    pc.add_argument("--cell-barcode-tag", default=None,
                    help="Per-cell breakdown via this BAM tag (e.g. CB). Opt-in; "
                         "no-ops with a warning if the tag is absent.")
    pc.add_argument("--per-gene-table", action="store_true",
                    help="Also write {sample}_per_gene.tsv with per-gene counts.")
    _add_common_indexing_args(pc)

    # basecomp
    pb = sub.add_parser("basecomp", help="Downstream base composition analysis")
    pb.add_argument("--bam", required=True)
    pb.add_argument("--gtf", required=True)
    pb.add_argument("--polya-sites", required=True)
    pb.add_argument("--fasta", required=True,
                    help="Indexed FASTA (gz OK with .fai + .gzi)")
    pb.add_argument("--sample", default="sample")
    pb.add_argument("--out-dir", default=".")
    pb.add_argument("--window", type=int, default=20,
                    help="Window size (nt) downstream of cleavage (default 20)")
    _add_common_indexing_args(pb)

    # compare
    pcmp = sub.add_parser("compare", help="Cross-sample comparison plots")
    pcmp.add_argument("--mode", choices=("classify", "basecomp"), required=True)
    pcmp.add_argument("--inputs", required=True,
                      help="Comma-separated JSON files to compare.")
    pcmp.add_argument("--out-dir", default=".")
    pcmp.add_argument("--verbose", action="store_true")

    # download-atlas
    pdl = sub.add_parser("download-atlas",
                         help="Download PolyASite 3.0 (and optionally GENCODE GTF)")
    pdl.add_argument("--genome", choices=("GRCh38", "GRCm39"), required=True)
    pdl.add_argument("--out-dir", default=".")
    pdl.add_argument("--url", default=None,
                     help="Override polyA atlas URL (otherwise built-in).")
    pdl.add_argument("--expected-sha256", default=None,
                     help="Pin a known SHA256 to verify reproducibility.")
    pdl.add_argument("--gtf-version", type=int, default=None,
                     help="If set, also fetch matching GENCODE GTF release "
                          "(human integer e.g. 45, mouse integer e.g. 35).")
    pdl.add_argument("--verbose", action="store_true")

    return parser


def _cmd_classify(args):
    from tecap.bam import analyse_bam, per_cell_summary, write_per_gene_table
    from tecap.gtf import build_gene_index
    from tecap.io import write_json
    from tecap.multiqc import build_mqc_payload, write_mqc_json
    from tecap.plotting import plot_single
    from tecap.polya import build_polya_index

    log = logging.getLogger("tecap.classify")
    log.info("sample=%s bam=%s", args.sample, args.bam)
    os.makedirs(args.out_dir, exist_ok=True)

    gene_index, gene_records = build_gene_index(args.gtf)
    allowed_types = (None if args.polya_types.lower() == "all"
                     else set(t.strip() for t in args.polya_types.split(",") if t.strip()))
    polya_index = build_polya_index(args.polya_sites,
                                    slop=args.polya_slop,
                                    allowed_types=allowed_types)

    results = analyse_bam(
        bam_path      = args.bam,
        gene_index    = gene_index,
        gene_records  = gene_records,
        polya_index   = polya_index,
        cov_threshold = args.coverage_threshold,
        min_mapq      = args.min_mapq,
        threads       = args.threads,
        cb_tag        = args.cell_barcode_tag,
    )
    results["sample"] = args.sample

    pas = results["pas_counts"]
    orient = results["orient"]
    orient_total = orient["match"] + orient["mismatch"]
    summary = {
        "sample":               args.sample,
        "platform":             args.platform,
        "total":                results["total"],
        "counts":               results["counts"],
        "fractions":            results["fractions"],
        "polya_slop":           args.polya_slop,
        "polya_types":          args.polya_types,
        "pas_counts":           pas,
        "pas_fraction":         {k: (v["pas"] / (v["pas"] + v["nopas"])
                                      if (v["pas"] + v["nopas"]) else None)
                                 for k, v in pas.items()},
        "orient_match_frac":    (orient["match"] / orient_total
                                 if orient_total else None),
        "utr_bin_labels":       UTR_BIN_LABELS,
        "utr_bin_total":        results["utr_bin_total"],
        "utr_bin_mecha_correct": results["utr_bin_mecha_correct"],
        "utr_bin_captured":     results["utr_bin_captured"],
        "utr_bin_mecha_rate":   results["utr_bin_mecha_rate"],
        "utr_bin_capt_rate":    results["utr_bin_capt_rate"],
        "cov_median_by_cat":    {k: float(np.median(v))
                                 for k, v in results["cov_fracs"].items() if v},
        "read_len_median_mecha": float(np.median(results["mecha_read_lengths"]))
                                  if results["mecha_read_lengths"] else None,
        "read_len_median_capt":  float(np.median(results["capt_read_lengths"]))
                                  if results["capt_read_lengths"] else None,
    }

    if args.cell_barcode_tag is not None:
        summary["cell_barcode_tag"] = args.cell_barcode_tag
        summary["per_cell_metrics"] = per_cell_summary(results["per_cell"])

    json_path = os.path.join(args.out_dir, f"{args.sample}_terminal_exon.json")

    if args.per_gene_table:
        tsv_path = os.path.join(args.out_dir, f"{args.sample}_per_gene.tsv")
        write_per_gene_table(results["per_gene"], gene_records, tsv_path)
        summary["per_gene_table_path"] = os.path.basename(tsv_path)
        log.info("Per-gene TSV: %s", tsv_path)

    write_json(json_path, summary)
    log.info("JSON: %s", json_path)

    mqc_path = os.path.join(args.out_dir, f"{args.sample}_tecap_mqc.json")
    write_mqc_json(mqc_path, build_mqc_payload(args.sample, results, summary))
    log.info("MultiQC: %s", mqc_path)

    plot_single(results, args.sample, args.out_dir, args.coverage_threshold)

    total = results["total"]
    print()
    print(f"  {'Category':<52} {'Count':>9}  {'%':>5}")
    print("  " + "-" * 72)
    for cat in CATEGORIES:
        c = results["counts"].get(cat, 0)
        pct = c / total * 100 if total > 0 else 0
        print(f"  {cat:<52} {c:>9,}  {pct:>4.1f}%")
    print("  " + "-" * 72)
    print(f"  {'TOTAL':<52} {total:>9,}  100.0%")
    print()
    print(f"  Median read length -- Captured      : {summary['read_len_median_capt']} bp")
    print(f"  Median read length -- MechA-correct : {summary['read_len_median_mecha']} bp")
    if summary["orient_match_frac"] is not None:
        print(f"  BAM/gene-strand match fraction      : {summary['orient_match_frac']:.4f}")
    print()
    print("  PAS hexamer split (PAS+ = atlas cluster has canonical AAUAAA-like):")
    for cat in (MECH_A_CORRECT, MECH_B_APA):
        p = pas[cat]["pas"]
        np_ = pas[cat]["nopas"]
        tot = p + np_
        if tot:
            print(f"    {cat:<55} PAS+ {p:>9,} ({p/tot*100:5.1f}%)  "
                  f"PAS- {np_:>9,} ({np_/tot*100:5.1f}%)")


def _cmd_basecomp(args):
    from tecap.basecomp import analyse, summarise
    from tecap.constants import BUCKETS
    from tecap.gtf import build_gene_index
    from tecap.io import write_json
    from tecap.plotting import plot_histograms_basecomp
    from tecap.polya import build_polya_index

    log = logging.getLogger("tecap.basecomp")
    log.info("sample=%s bam=%s fasta=%s window=%d",
             args.sample, args.bam, args.fasta, args.window)
    os.makedirs(args.out_dir, exist_ok=True)

    gene_index, gene_records = build_gene_index(args.gtf)
    allowed_types = (None if args.polya_types.lower() == "all"
                     else set(t.strip() for t in args.polya_types.split(",") if t.strip()))
    polya_index = build_polya_index(args.polya_sites,
                                    slop=args.polya_slop,
                                    allowed_types=allowed_types)

    pct_a_by_bucket, _ = analyse(
        bam_path      = args.bam,
        fasta_path    = args.fasta,
        gene_index    = gene_index,
        gene_records  = gene_records,
        polya_index   = polya_index,
        cov_threshold = args.coverage_threshold,
        min_mapq      = args.min_mapq,
        window        = args.window,
        threads       = args.threads,
    )

    summary = summarise(pct_a_by_bucket)
    summary["sample"]      = args.sample
    summary["platform"]    = args.platform
    summary["window"]      = args.window
    summary["polya_slop"]  = args.polya_slop
    summary["polya_types"] = args.polya_types

    json_path = os.path.join(args.out_dir, f"{args.sample}_basecomp.json")
    write_json(json_path, summary)
    log.info("JSON: %s", json_path)

    plot_histograms_basecomp(summary, args.sample, args.out_dir, args.window)

    print()
    print(f"  {'Bucket':<20} {'n':>12}  {'median %A':>10}  {'>=60%':>7}  {'30-50%':>7}")
    print("  " + "-" * 64)
    for b in BUCKETS:
        d = summary["buckets"][b]
        if d["n"] == 0:
            print(f"  {b:<20} {'0':>12}  {'-':>10}  {'-':>7}  {'-':>7}")
            continue
        print(f"  {b:<20} {d['n']:>12,}  {d['median']*100:>9.1f}%  "
              f"{d['frac_ge_60']*100:>6.1f}%  {d['frac_30_50']*100:>6.1f}%")


def _cmd_compare(args):
    from tecap.compare import run_compare

    os.makedirs(args.out_dir, exist_ok=True)
    inputs = [p.strip() for p in args.inputs.split(",") if p.strip()]
    run_compare(args.mode, inputs, args.out_dir)


def _cmd_download_atlas(args):
    from tecap.download import fetch_gencode_gtf, fetch_polya_atlas

    log = logging.getLogger("tecap.download")
    bed_path, bed_sha = fetch_polya_atlas(
        genome=args.genome, out_dir=args.out_dir,
        expected_sha256=args.expected_sha256, url=args.url,
    )
    log.info("PolyASite atlas: %s (sha256=%s)", bed_path, bed_sha)

    if args.gtf_version is not None:
        gtf_path, gtf_sha = fetch_gencode_gtf(
            genome=args.genome, version=args.gtf_version, out_dir=args.out_dir,
        )
        log.info("GENCODE GTF: %s (sha256=%s)", gtf_path, gtf_sha)


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)

    _setup_logging(getattr(args, "verbose", False))

    if args.cmd in ("classify", "basecomp"):
        _check_platform(args.platform)

    if args.cmd == "classify":
        _cmd_classify(args)
    elif args.cmd == "basecomp":
        _cmd_basecomp(args)
    elif args.cmd == "compare":
        _cmd_compare(args)
    elif args.cmd == "download-atlas":
        _cmd_download_atlas(args)
    else:
        parser.error(f"Unknown command: {args.cmd}")


if __name__ == "__main__":
    main()
