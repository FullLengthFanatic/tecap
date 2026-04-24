"""Base composition analysis in the N-nt reference window downstream of each
read's 3' cleavage site (mRNA orientation)."""

import logging
import multiprocessing as mp
from collections import defaultdict

import numpy as np
import pysam

from tecap.classify import classify_read
from tecap.constants import (
    BUCKETS, BUCKET_CAPT, BUCKET_MA_PAS, BUCKET_MA_NOPAS,
    BUCKET_MB_APA_PAS, BUCKET_MB_APA_NOPAS, BUCKET_MB_EXON,
    BUCKET_MB_ASPECI, BUCKET_IP_CDS,
    CAPTURED, MECH_A_CORRECT, MECH_B_APA, MECH_B_EXON, MECH_B_ASPECI,
    INTERNAL_PRIME_TE_CDS,
)

log = logging.getLogger(__name__)


_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


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


def pct_a(seq):
    if not seq:
        return None
    return seq.count("A") / len(seq)


def _process_chrom_bc(args):
    (bam_path, fasta_path, chrom, gene_index, gene_records, polya_index,
     cov_threshold, min_mapq, window) = args

    pct_a_by_bucket = defaultdict(list)
    counts          = defaultdict(int)
    n_no_seq        = 0

    try:
        bam   = pysam.AlignmentFile(bam_path, "rb")
        fasta = pysam.FastaFile(fasta_path)
    except Exception as e:
        log.error("Failed to open BAM or FASTA for contig %s: %s", chrom, e)
        return pct_a_by_bucket, counts, 0

    try:
        iterator = bam.fetch(contig=chrom)
    except ValueError:
        bam.close()
        fasta.close()
        return pct_a_by_bucket, counts, 0

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
        cat, _, r3p, has_pas = classify_read(read, rec, polya_index, cov_threshold)
        bucket = assign_bucket(cat, has_pas)
        if bucket is None:
            continue

        seq = downstream_window(fasta, chrom, r3p, rec["strand"], window)
        a = pct_a(seq)
        if a is None:
            n_no_seq += 1
            continue

        pct_a_by_bucket[bucket].append(a)
        counts[bucket] += 1

    bam.close()
    fasta.close()
    return pct_a_by_bucket, counts, n_no_seq


def analyse(bam_path, fasta_path, gene_index, gene_records, polya_index,
            cov_threshold, min_mapq, window, threads=1):
    bam = pysam.AlignmentFile(bam_path, "rb")
    contigs = list(bam.references)
    bam.close()

    log.info("Basecomp analysing %d contigs with %d worker(s) (window=%d)",
             len(contigs), threads, window)

    worker_args = [
        (bam_path, fasta_path, c, gene_index, gene_records, polya_index,
         cov_threshold, min_mapq, window)
        for c in contigs
    ]

    pct_a_by_bucket = defaultdict(list)
    bucket_counts   = defaultdict(int)
    n_no_seq_total  = 0

    if threads <= 1:
        results_iter = (_process_chrom_bc(wa) for wa in worker_args)
    else:
        ctx = mp.get_context("fork")
        pool = ctx.Pool(processes=threads)
        results_iter = pool.imap_unordered(_process_chrom_bc, worker_args)

    for pab, cnt, nns in results_iter:
        for b, vs in pab.items():
            pct_a_by_bucket[b].extend(vs)
        for b, v in cnt.items():
            bucket_counts[b] += v
        n_no_seq_total += nns

    if threads > 1:
        pool.close()
        pool.join()

    log.info("Basecomp done; %s reads dropped for missing/short window",
             f"{n_no_seq_total:,}")

    return pct_a_by_bucket, dict(bucket_counts)


def summarise(pct_a_by_bucket):
    hist_bins = np.linspace(0.0, 1.0, 21)  # 5% bins
    out = {"bin_edges": hist_bins.tolist(), "buckets": {}}
    for b in BUCKETS:
        vals = np.asarray(pct_a_by_bucket.get(b, []), dtype=np.float32)
        if len(vals) == 0:
            out["buckets"][b] = {
                "n": 0, "median": None, "mean": None,
                "frac_ge_60":  None, "frac_30_50": None, "frac_ge_40": None,
                "hist": [0] * (len(hist_bins) - 1),
            }
            continue
        hist, _ = np.histogram(vals, bins=hist_bins)
        out["buckets"][b] = {
            "n":           int(len(vals)),
            "median":      float(np.median(vals)),
            "mean":        float(np.mean(vals)),
            "frac_ge_60":  float((vals >= 0.60).mean()),
            "frac_30_50":  float(((vals >= 0.30) & (vals < 0.50)).mean()),
            "frac_ge_40":  float((vals >= 0.40).mean()),
            "hist":        hist.tolist(),
        }
    return out
