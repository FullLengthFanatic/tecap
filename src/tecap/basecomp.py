"""Base composition analysis in the N-nt reference window downstream of each
read's 3' cleavage site (mRNA orientation)."""

import logging
import multiprocessing as mp
from collections import defaultdict

import numpy as np
import pysam

from tecap.bam import choose_gene
from tecap.classify import classify_read
from tecap.constants import (
    BASECOMP_GE_HI,
    BASECOMP_GE_MID,
    BASECOMP_HIST_N,
    BASECOMP_MID_HI,
    BASECOMP_MID_LO,
    BUCKET_CAPT,
    BUCKET_IP_CDS,
    BUCKET_MA_NOPAS,
    BUCKET_MA_PAS,
    BUCKET_MB_APA_NOPAS,
    BUCKET_MB_APA_PAS,
    BUCKET_MB_ASPECI,
    BUCKET_MB_EXON,
    BUCKETS,
    CAPTURED,
    DEFINITION_SOURCE,
    INTERNAL_PRIME_TE_CDS,
    LABEL_WINDOW,
    LABELS,
    MECH_A_CORRECT,
    MECH_B_APA,
    MECH_B_ASPECI,
    MECH_B_EXON,
    NONTEMPLATED_MAX,
    SITE_ANCHOR_DEFAULT,
    SITE_RUN_MIN,
    SITE_RUN_TOLERANCE,
    TEMPLATED_MIN,
)
from tecap.site import label_of, read_site

log = logging.getLogger(__name__)


_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


_WORKER = {
    "bam_path":      None,
    "fasta_path":    None,
    "gene_index":    None,
    "gene_records":  None,
    "polya_index":   None,
    "cov_threshold": None,
    "min_mapq":      None,
    "window":        None,
    "site_anchor":   None,
    "skip_antisense": None,
}


def rc(seq):
    return seq.translate(_COMP)[::-1]


def assign_bucket(cat, has_pas):
    """Map classify_read output to a base-composition bucket. Returns None for
    buckets we don't analyse (MECH_A_INTERNAL_UTR, MECH_A_NO_CDS, MECH_C)."""
    if cat == CAPTURED:
        return BUCKET_CAPT
    if cat == MECH_A_CORRECT:
        return BUCKET_MA_PAS if has_pas else BUCKET_MA_NOPAS
    if cat == MECH_B_APA:
        return BUCKET_MB_APA_PAS if has_pas else BUCKET_MB_APA_NOPAS
    if cat == MECH_B_EXON:
        return BUCKET_MB_EXON
    if cat == MECH_B_ASPECI:
        return BUCKET_MB_ASPECI
    if cat == INTERNAL_PRIME_TE_CDS:
        return BUCKET_IP_CDS
    return None


def downstream_window(fasta, chrom, pos, strand, win):
    """Return the reference sequence of length `win` immediately 3' of `pos` on
    the mRNA strand. For + genes: [pos, pos+win). For -: RC of [pos-win, pos).
    Returns '' on failure or out-of-bounds."""
    try:
        if strand == "+":
            seq = fasta.fetch(chrom, pos, pos + win)
        else:
            start = pos - win
            if start < 0:
                return ""
            seq = rc(fasta.fetch(chrom, start, pos))
    except (KeyError, ValueError):
        return ""
    seq = seq.upper()
    if len(seq) != win:
        return ""
    return seq


def _process_chrom_bc(chrom):
    bam_path      = _WORKER["bam_path"]
    fasta_path    = _WORKER["fasta_path"]
    gene_index    = _WORKER["gene_index"]
    gene_records  = _WORKER["gene_records"]
    polya_index   = _WORKER["polya_index"]
    cov_threshold = _WORKER["cov_threshold"]
    min_mapq      = _WORKER["min_mapq"]
    window        = _WORKER["window"]
    site_anchor   = _WORKER["site_anchor"]
    skip_antisense = _WORKER["skip_antisense"]

    a_count_by_bucket = defaultdict(list)
    counts            = defaultdict(int)
    n_no_seq          = 0
    orient            = {"match": 0, "mismatch": 0}
    site_stats        = {"observed": 0, "fallback": 0, "antisense_skipped": 0,
                         "antisense_reassigned": 0, "labels": {}}

    try:
        bam   = pysam.AlignmentFile(bam_path, "rb")
        fasta = pysam.FastaFile(fasta_path)
    except Exception as e:
        log.error("Failed to open BAM or FASTA for contig %s: %s", chrom, e)
        return a_count_by_bucket, counts, 0, orient, site_stats

    try:
        iterator = bam.fetch(contig=chrom)
    except ValueError:
        bam.close()
        fasta.close()
        return a_count_by_bucket, counts, 0, orient, site_stats

    genes = gene_index[chrom]
    for read in iterator:
        best_gid, status = choose_gene(read, genes, gene_records, min_mapq,
                                       skip_antisense)
        if status is None:
            continue

        orient["match" if status == "sense" else "mismatch"] += 1
        if status == "skipped":
            site_stats["antisense_skipped"] += 1
            continue
        if status == "reassigned":
            site_stats["antisense_reassigned"] += 1

        rec = gene_records[best_gid]
        site, observed = read_site(read, rec["strand"], site_anchor)
        site_stats["observed" if observed else "fallback"] += 1
        cat, _, _, has_pas = classify_read(read, rec, polya_index, cov_threshold,
                                           read_3prime=site)
        bucket = assign_bucket(cat, has_pas)
        if bucket is None:
            continue

        seq = downstream_window(fasta, chrom, site, rec["strand"], window)
        if not seq:
            n_no_seq += 1
            continue

        k = seq.count("A")
        a_count_by_bucket[bucket].append(k)
        counts[bucket] += 1
        if observed and window == LABEL_WINDOW:
            by_label = site_stats["labels"].setdefault(bucket, dict.fromkeys(LABELS, 0))
            by_label[label_of(k)] += 1

    bam.close()
    fasta.close()
    return a_count_by_bucket, counts, n_no_seq, orient, site_stats


def _populated_contigs(bam_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        try:
            stats = bam.get_index_statistics()
            return [s.contig for s in stats if s.mapped > 0]
        except (AttributeError, ValueError):
            return list(bam.references)
    finally:
        bam.close()


def analyse(bam_path, fasta_path, gene_index, gene_records, polya_index,
            cov_threshold, min_mapq, window, threads=1,
            site_anchor=SITE_ANCHOR_DEFAULT, skip_antisense=False):
    """Returns (a_count_by_bucket, bucket_counts, orient, site_stats).

    a_count_by_bucket holds, per bucket, the integer adenine count of each
    read's downstream window. site_anchor picks the read 3' end used for
    classification and as the start of the downstream window (see
    tecap.site.read_site); site_stats counts the reads whose site was
    observed and those that fell back to the alignment end, and under
    "labels" the label of each bucketed read with an observed site and a
    window, per bucket. Labels are counted only when window is LABEL_WINDOW.

    orient compares each read's strand with its best-overlap gene on either
    strand. With skip_antisense (stranded input), a read antisense to that
    gene is bucketed against the best-overlap gene on its own strand and
    counted in site_stats["antisense_reassigned"], or, when no gene on its
    strand overlaps it, counted in orient and in
    site_stats["antisense_skipped"] and gets no site, bucket or window (see
    tecap.bam.choose_gene).
    """
    contigs = _populated_contigs(bam_path)
    log.info("Basecomp analysing %d populated contigs with %d worker(s) (window=%d)",
             len(contigs), threads, window)

    _WORKER["bam_path"]      = bam_path
    _WORKER["fasta_path"]    = fasta_path
    _WORKER["gene_index"]    = gene_index
    _WORKER["gene_records"]  = gene_records
    _WORKER["polya_index"]   = polya_index
    _WORKER["cov_threshold"] = cov_threshold
    _WORKER["min_mapq"]      = min_mapq
    _WORKER["window"]        = window
    _WORKER["site_anchor"]   = site_anchor
    _WORKER["skip_antisense"] = skip_antisense

    a_count_by_bucket = defaultdict(list)
    bucket_counts     = defaultdict(int)
    n_no_seq_total    = 0
    orient            = {"match": 0, "mismatch": 0}
    site_stats        = {"observed": 0, "fallback": 0, "antisense_skipped": 0,
                         "antisense_reassigned": 0, "labels": {}}

    try:
        if threads <= 1:
            results_iter = (_process_chrom_bc(c) for c in contigs)
            for acb, cnt, nns, ori, sst in results_iter:
                for b, vs in acb.items():
                    a_count_by_bucket[b].extend(vs)
                for b, v in cnt.items():
                    bucket_counts[b] += v
                n_no_seq_total += nns
                orient["match"]    += ori["match"]
                orient["mismatch"] += ori["mismatch"]
                site_stats["observed"] += sst["observed"]
                site_stats["fallback"] += sst["fallback"]
                site_stats["antisense_skipped"] += sst["antisense_skipped"]
                site_stats["antisense_reassigned"] += sst["antisense_reassigned"]
                for b, by_label in sst["labels"].items():
                    tot = site_stats["labels"].setdefault(b, dict.fromkeys(LABELS, 0))
                    for x, v in by_label.items():
                        tot[x] += v
        else:
            ctx = mp.get_context("fork")
            with ctx.Pool(processes=threads) as pool:
                for acb, cnt, nns, ori, sst in pool.imap_unordered(_process_chrom_bc, contigs):
                    for b, vs in acb.items():
                        a_count_by_bucket[b].extend(vs)
                    for b, v in cnt.items():
                        bucket_counts[b] += v
                    n_no_seq_total += nns
                    orient["match"]    += ori["match"]
                    orient["mismatch"] += ori["mismatch"]
                    site_stats["observed"] += sst["observed"]
                    site_stats["fallback"] += sst["fallback"]
                    site_stats["antisense_skipped"] += sst["antisense_skipped"]
                    site_stats["antisense_reassigned"] += sst["antisense_reassigned"]
                    for b, by_label in sst["labels"].items():
                        tot = site_stats["labels"].setdefault(b, dict.fromkeys(LABELS, 0))
                        for x, v in by_label.items():
                            tot[x] += v
    finally:
        for k in _WORKER:
            _WORKER[k] = None

    log.info("Basecomp done; %s reads dropped for missing/short window",
             f"{n_no_seq_total:,}")

    return a_count_by_bucket, dict(bucket_counts), orient, site_stats


def summarise(a_count_by_bucket, window=LABEL_WINDOW, labels=None):
    """Per-bucket summary of the adenine counts in the `window` nt window.

    The adenine share bands (median, mean, hist, frac_*) are computed on
    count / window. The label fields come from `labels` (site_stats["labels"]
    of analyse); they are None when labels is None or window is not
    LABEL_WINDOW, and frac_<label> is None when the bucket has no read with
    a site.
    """
    hist_bins = np.linspace(0.0, 1.0, BASECOMP_HIST_N)  # 5% bins
    labelled = labels is not None and window == LABEL_WINDOW
    out = {
        "definition_source": DEFINITION_SOURCE,
        "bin_edges": hist_bins.tolist(),
        "thresholds": {
            "frac_ge_60": BASECOMP_GE_HI,
            "frac_30_50": [BASECOMP_MID_LO, BASECOMP_MID_HI],
            "frac_ge_40": BASECOMP_GE_MID,
            "hist_bins":  BASECOMP_HIST_N,
            "templated_min":      TEMPLATED_MIN,
            "nontemplated_max":   NONTEMPLATED_MAX,
            "label_window":       LABEL_WINDOW,
            "site_run_min":       SITE_RUN_MIN,
            "site_run_tolerance": SITE_RUN_TOLERANCE,
        },
        "buckets": {},
    }
    for b in BUCKETS:
        counts = np.asarray(a_count_by_bucket.get(b, []), dtype=np.float64)
        vals = (counts / window).astype(np.float32)
        if labelled:
            n = {x: labels.get(b, {}).get(x, 0) for x in LABELS}
            n_site = sum(n.values())
            label_fields = {
                "n_site": n_site,
                **{f"n_{x}": n[x] for x in LABELS},
                **{f"frac_{x}": n[x] / n_site if n_site else None for x in LABELS},
            }
        else:
            label_fields = dict.fromkeys(
                ["n_site", *(f"n_{x}" for x in LABELS), *(f"frac_{x}" for x in LABELS)])
        if len(vals) == 0:
            out["buckets"][b] = {
                "n": 0, "median": None, "mean": None,
                "frac_ge_60":  None, "frac_30_50": None, "frac_ge_40": None,
                "hist": [0] * (len(hist_bins) - 1),
                **label_fields,
            }
            continue
        hist, _ = np.histogram(vals, bins=hist_bins)
        out["buckets"][b] = {
            "n":           int(len(vals)),
            "median":      float(np.median(vals)),
            "mean":        float(np.mean(vals)),
            "frac_ge_60":  float((vals >= BASECOMP_GE_HI).mean()),
            "frac_30_50":  float(((vals >= BASECOMP_MID_LO) & (vals < BASECOMP_MID_HI)).mean()),
            "frac_ge_40":  float((vals >= BASECOMP_GE_MID).mean()),
            "hist":        hist.tolist(),
            **label_fields,
        }
    return out
