"""BAM iteration, classification dispatch, accumulator merging.

Chromosomes are the parallelism unit. Gene and polyA indices are placed in
module-level globals before the worker pool is forked, so each worker inherits
them via copy-on-write rather than receiving them through pickle. Contigs with
no mapped reads (per the BAM index) are skipped at the parent."""

import csv
import logging
import multiprocessing as mp
from collections import defaultdict

import pysam

from tecap.classify import classify_read
from tecap.constants import (
    CAPTURED,
    CATEGORIES,
    MECH_A_CORRECT,
    MECH_B_APA,
    UTR_BIN_LABELS,
)
from tecap.polya import utr_bin

log = logging.getLogger(__name__)


_WORKER = {
    "bam_path":      None,
    "gene_index":    None,
    "gene_records":  None,
    "polya_index":   None,
    "cov_threshold": None,
    "min_mapq":      None,
    "cb_tag":        None,
}


def _empty_accumulator():
    n_bins = len(UTR_BIN_LABELS)
    return {
        "counts":               defaultdict(int),
        "cov_fracs":             defaultdict(list),
        "pas_counts":           {MECH_A_CORRECT: {"pas": 0, "nopas": 0},
                                  MECH_B_APA:     {"pas": 0, "nopas": 0}},
        "orient":               {"match": 0, "mismatch": 0},
        "mecha_read_lengths":    [],
        "capt_read_lengths":     [],
        "mecha_cov_scatter":     [],
        "utr_bin_total":         [0] * n_bins,
        "utr_bin_mecha_correct": [0] * n_bins,
        "utr_bin_captured":      [0] * n_bins,
        "per_cell":              {},
        "per_gene":              {},
        "cb_seen":               0,
        "cb_missing":            0,
    }


def _bump(d, k1, k2):
    """Bump d[k1][k2] by 1, creating the inner dict on demand. Keeps nested
    counters picklable (plain dicts, no defaultdict factories)."""
    inner = d.get(k1)
    if inner is None:
        inner = {}
        d[k1] = inner
    inner[k2] = inner.get(k2, 0) + 1


def _merge_accumulator(dst, src):
    for k, v in src["counts"].items():
        dst["counts"][k] += v
    for k, vs in src["cov_fracs"].items():
        dst["cov_fracs"][k].extend(vs)
    for cat in dst["pas_counts"]:
        dst["pas_counts"][cat]["pas"]   += src["pas_counts"][cat]["pas"]
        dst["pas_counts"][cat]["nopas"] += src["pas_counts"][cat]["nopas"]
    dst["orient"]["match"]    += src["orient"]["match"]
    dst["orient"]["mismatch"] += src["orient"]["mismatch"]
    dst["mecha_read_lengths"].extend(src["mecha_read_lengths"])
    dst["capt_read_lengths"].extend(src["capt_read_lengths"])
    dst["mecha_cov_scatter"].extend(src["mecha_cov_scatter"])
    for i in range(len(dst["utr_bin_total"])):
        dst["utr_bin_total"][i]         += src["utr_bin_total"][i]
        dst["utr_bin_mecha_correct"][i] += src["utr_bin_mecha_correct"][i]
        dst["utr_bin_captured"][i]      += src["utr_bin_captured"][i]
    for cb, d in src["per_cell"].items():
        dst_inner = dst["per_cell"].setdefault(cb, {})
        for k, v in d.items():
            dst_inner[k] = dst_inner.get(k, 0) + v
    for gid, d in src["per_gene"].items():
        dst_inner = dst["per_gene"].setdefault(gid, {})
        for k, v in d.items():
            dst_inner[k] = dst_inner.get(k, 0) + v
    dst["cb_seen"]    += src["cb_seen"]
    dst["cb_missing"] += src["cb_missing"]


def _process_chrom(chrom):
    """Worker: process all reads in one contig, return an accumulator. Reads
    its inputs from the module-level _WORKER dict (set by the parent before
    forking the pool)."""
    bam_path      = _WORKER["bam_path"]
    gene_index    = _WORKER["gene_index"]
    gene_records  = _WORKER["gene_records"]
    polya_index   = _WORKER["polya_index"]
    cov_threshold = _WORKER["cov_threshold"]
    min_mapq      = _WORKER["min_mapq"]
    cb_tag        = _WORKER["cb_tag"]

    acc = _empty_accumulator()
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        log.error("Failed to open BAM for contig %s: %s", chrom, e)
        return acc

    try:
        iterator = bam.fetch(contig=chrom)
    except ValueError:
        bam.close()
        return acc

    for read in iterator:
        if (read.is_unmapped or read.is_secondary or
                read.is_supplementary or not read.cigartuples):
            continue
        if read.mapping_quality < min_mapq:
            continue
        if not any(op == 3 for op, _ in read.cigartuples):
            continue

        hits = gene_index[chrom][read.reference_start:read.reference_end]
        if not hits:
            continue

        best_gid, best_ovlp = None, 0
        for iv in hits:
            gid = iv.data
            rec = gene_records[gid]
            if rec["strand"] not in ("+", "-"):
                continue
            ovlp = (min(read.reference_end, iv.end) -
                    max(read.reference_start, iv.begin))
            if ovlp > best_ovlp:
                best_ovlp, best_gid = ovlp, gid

        if best_gid is None:
            continue

        rec = gene_records[best_gid]
        cat, cov, _, has_pas = classify_read(read, rec, polya_index, cov_threshold)
        acc["counts"][cat] += 1
        acc["cov_fracs"][cat].append(cov)

        if cat in acc["pas_counts"]:
            acc["pas_counts"][cat]["pas" if has_pas else "nopas"] += 1

        expect_reverse = (rec["strand"] == "-")
        if read.is_reverse == expect_reverse:
            acc["orient"]["match"] += 1
        else:
            acc["orient"]["mismatch"] += 1

        ub = utr_bin(rec["utr_length"])
        acc["utr_bin_total"][ub] += 1
        if cat == MECH_A_CORRECT:
            acc["utr_bin_mecha_correct"][ub] += 1
        elif cat == CAPTURED:
            acc["utr_bin_captured"][ub] += 1

        rlen = read.query_alignment_length
        if cat == MECH_A_CORRECT:
            acc["mecha_read_lengths"].append(rlen)
            acc["mecha_cov_scatter"].append((rlen, cov))
        elif cat == CAPTURED:
            acc["capt_read_lengths"].append(rlen)

        _bump(acc["per_gene"], best_gid, cat)

        if cb_tag is not None:
            try:
                cb = read.get_tag(cb_tag)
                _bump(acc["per_cell"], cb, cat)
                acc["cb_seen"] += 1
            except KeyError:
                acc["cb_missing"] += 1

    bam.close()
    return acc


def _populated_contigs(bam_path):
    """Return contigs with mapped > 0 according to the BAM index. Falls back
    to the full reference list if the index lacks per-contig stats."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        try:
            stats = bam.get_index_statistics()
            return [s.contig for s in stats if s.mapped > 0]
        except (AttributeError, ValueError):
            return list(bam.references)
    finally:
        bam.close()


def analyse_bam(bam_path, gene_index, gene_records, polya_index,
                cov_threshold, min_mapq, threads=1, cb_tag=None,
                cb_probe_reads=10_000):
    """Process a BAM in parallel by contig. Returns a results dict matching
    the v0.1 JSON schema.

    Empty contigs (mapped=0 in the BAM index) are skipped. With threads > 1,
    the pool is forked from a parent that has placed the gene/polyA indices
    in module-level globals, so workers inherit them via copy-on-write rather
    than receive them through pickle.
    """

    probed_cb_tag = cb_tag
    if cb_tag is not None:
        probed_cb_tag = _probe_cb_tag(bam_path, cb_tag, cb_probe_reads)

    contigs = _populated_contigs(bam_path)
    log.info("Analysing %d populated contigs with %d worker(s)",
             len(contigs), threads)

    _WORKER["bam_path"]      = bam_path
    _WORKER["gene_index"]    = gene_index
    _WORKER["gene_records"]  = gene_records
    _WORKER["polya_index"]   = polya_index
    _WORKER["cov_threshold"] = cov_threshold
    _WORKER["min_mapq"]      = min_mapq
    _WORKER["cb_tag"]        = probed_cb_tag

    merged = _empty_accumulator()

    try:
        if threads <= 1:
            for chrom in contigs:
                acc = _process_chrom(chrom)
                _merge_accumulator(merged, acc)
        else:
            ctx = mp.get_context("fork")
            with ctx.Pool(processes=threads) as pool:
                for acc in pool.imap_unordered(_process_chrom, contigs):
                    _merge_accumulator(merged, acc)
    finally:
        for k in _WORKER:
            _WORKER[k] = None

    if cb_tag is not None and probed_cb_tag is None:
        log.warning("Cell barcode tag %r not found in first %d reads; "
                    "per-cell output will be empty", cb_tag, cb_probe_reads)

    total     = sum(merged["counts"].values())
    fractions = {k: v / total if total > 0 else 0 for k, v in merged["counts"].items()}

    n_bins = len(UTR_BIN_LABELS)
    utr_bin_mecha_rate = [
        merged["utr_bin_mecha_correct"][i] / merged["utr_bin_total"][i]
        if merged["utr_bin_total"][i] > 0 else 0
        for i in range(n_bins)
    ]
    utr_bin_capt_rate = [
        merged["utr_bin_captured"][i] / merged["utr_bin_total"][i]
        if merged["utr_bin_total"][i] > 0 else 0
        for i in range(n_bins)
    ]

    log.info("Processed %s classified reads", f"{total:,}")

    return {
        "counts":                dict(merged["counts"]),
        "fractions":             fractions,
        "total":                 total,
        "cov_fracs":             dict(merged["cov_fracs"]),
        "pas_counts":            merged["pas_counts"],
        "orient":                merged["orient"],
        "mecha_read_lengths":    merged["mecha_read_lengths"],
        "capt_read_lengths":     merged["capt_read_lengths"],
        "mecha_cov_scatter":     merged["mecha_cov_scatter"],
        "utr_bin_total":         merged["utr_bin_total"],
        "utr_bin_mecha_correct": merged["utr_bin_mecha_correct"],
        "utr_bin_captured":      merged["utr_bin_captured"],
        "utr_bin_mecha_rate":    utr_bin_mecha_rate,
        "utr_bin_capt_rate":     utr_bin_capt_rate,
        "per_cell":              merged["per_cell"],
        "per_gene":              merged["per_gene"],
        "cb_tag_seen":           probed_cb_tag is not None,
    }


def _probe_cb_tag(bam_path, cb_tag, n_probe):
    """Return cb_tag if any of the first n_probe reads carry it, else None."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        found = False
        for i, read in enumerate(bam.fetch(until_eof=True)):
            if i >= n_probe:
                break
            try:
                read.get_tag(cb_tag)
                found = True
                break
            except KeyError:
                continue
        return cb_tag if found else None
    finally:
        bam.close()


def write_per_gene_table(per_gene, gene_records, out_path):
    """Write a TSV: gene_id, te_length, utr_length, and one count column per
    bucket in CATEGORIES."""
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        header = ["gene_id", "chrom", "strand", "te_length", "utr_length"] + CATEGORIES
        w.writerow(header)
        for gid, counts in per_gene.items():
            rec = gene_records.get(gid)
            if rec is None:
                continue
            row = [gid, rec["chrom"], rec["strand"],
                   rec["te_length"], rec["utr_length"]]
            row.extend(counts.get(cat, 0) for cat in CATEGORIES)
            w.writerow(row)


def per_cell_summary(per_cell):
    """Compact per-cell dict: {cb: {counts: {cat: n}, total: N, fractions: ...}}."""
    out = {}
    for cb, counts in per_cell.items():
        total = sum(counts.values())
        if total == 0:
            continue
        out[cb] = {
            "total":     total,
            "counts":    dict(counts),
            "fractions": {k: v / total for k, v in counts.items()},
        }
    return out
