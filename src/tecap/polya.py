"""Strand-aware polyA site atlas handling (PolyASite 3.0)."""

import logging
from collections import defaultdict

from intervaltree import IntervalTree

from tecap.constants import UTR_BINS
from tecap.gtf import _open_any

log = logging.getLogger(__name__)


def build_polya_index(bed_path, slop, allowed_types=None):
    """Strand-aware polyA index, filtered to allowed cluster types.

    Expects PolyASite 3.0 BED: col6=strand, col10=type, col11=PAS hexamer.
    Each interval carries a boolean has_pas so classify_read can tell PAS+
    clusters (real APA) from PAS- clusters (likely internal-priming artifacts
    in the atlas build). Returns dict[(chrom, strand)] -> IntervalTree[has_pas].
    """
    keep = ",".join(sorted(allowed_types)) if allowed_types else "ALL"
    log.info("Loading polyA atlas (strand-aware, +-%d bp, types=%s): %s",
             slop, keep, bed_path)
    index = defaultdict(IntervalTree)
    n_kept = n_wrong_type = n_no_strand = 0
    with _open_any(bed_path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 6:
                continue
            chrom  = parts[0]
            strand = parts[5]
            if strand not in ("+", "-"):
                n_no_strand += 1
                continue
            if allowed_types and len(parts) >= 10:
                ctype = parts[9]
                if ctype not in allowed_types:
                    n_wrong_type += 1
                    continue
            pas_field = parts[10] if len(parts) >= 11 else "NaN"
            has_pas   = pas_field not in ("NaN", "", ".", "nan")
            start = max(0, int(parts[1]) - slop)
            end   = int(parts[2]) + slop
            if end > start:
                index[(chrom, strand)][start:end] = has_pas
                n_kept += 1
    log.info("polyA atlas: %s kept, %s wrong-type skipped, %s no-strand skipped",
             f"{n_kept:,}", f"{n_wrong_type:,}", f"{n_no_strand:,}")
    return index


def hits_polya(chrom, pos, strand, polya_index):
    """Return (hit, has_pas). has_pas=True iff ANY overlapping cluster carries
    a canonical PAS hexamer."""
    ivs = polya_index[(chrom, strand)][pos:pos + 1]
    if not ivs:
        return False, False
    return True, any(iv.data for iv in ivs)


def in_upstream_exon(pos, gene_rec):
    return len(gene_rec["upstream_exon_tree"][pos:pos + 1]) > 0


def utr_bin(utr_len):
    """Return the UTR length bin index for a given UTR length."""
    for i in range(len(UTR_BINS) - 1):
        if UTR_BINS[i] <= utr_len < UTR_BINS[i + 1]:
            return i
    return len(UTR_BINS) - 2
