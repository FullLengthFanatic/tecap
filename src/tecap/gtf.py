"""GTF parsing: build a per-gene index with terminal-exon and UTR coordinates."""

import gzip
import logging
import re
from collections import defaultdict

from intervaltree import IntervalTree

log = logging.getLogger(__name__)


def _open_any(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def build_gene_index(gtf_path):
    """Parse a GTF and return (gene_index, gene_records).

    gene_index: dict[chrom] -> IntervalTree with gene_id as data.
    gene_records: dict[gene_id] -> record with te_start, te_end, te_length,
    cds_end_genomic, utr_length, upstream_exon_tree, etc.
    """
    log.info("Parsing GTF: %s", gtf_path)

    def _attr(field, s):
        m = re.search(rf'{field}\s+"([^"]+)"', s)
        return m.group(1) if m else None

    genes = {}
    with _open_any(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature not in ("exon", "CDS"):
                continue
            chrom  = parts[0]
            start  = int(parts[3]) - 1
            end    = int(parts[4])
            strand = parts[6]
            gid    = _attr("gene_id", parts[8])
            if gid is None:
                continue
            if gid not in genes:
                genes[gid] = {"chrom": chrom, "strand": strand,
                              "exons": [], "cds": []}
            if feature == "exon":
                genes[gid]["exons"].append((start, end))
            else:
                genes[gid]["cds"].append((start, end))

    log.info("%s genes parsed", f"{len(genes):,}")

    gene_index   = defaultdict(IntervalTree)
    gene_records = {}

    for gid, g in genes.items():
        if not g["exons"]:
            continue
        strand = g["strand"]
        chrom  = g["chrom"]
        exons  = sorted(set(g["exons"]))

        if strand == "+":
            max_end  = max(e for _, e in exons)
            te_exons = [e for e in exons if e[1] == max_end]
            te_start = min(s for s, _ in te_exons)
            te_end   = max_end
        else:
            min_start = min(s for s, _ in exons)
            te_exons  = [e for e in exons if e[0] == min_start]
            te_start  = min_start
            te_end    = max(e for _, e in te_exons)

        te_length = te_end - te_start

        cds_end_genomic = None
        if g["cds"]:
            if strand == "+":
                cds_end_genomic = max(e for _, e in g["cds"])
            else:
                cds_end_genomic = min(s for s, _ in g["cds"])

        gene_start = min(s for s, _ in exons)
        gene_end   = max(e for _, e in exons)

        upstream_exon_tree = IntervalTree()
        for es, ee in exons:
            if strand == "+" and ee == te_end:
                continue
            if strand == "-" and es == te_start:
                continue
            if ee > es:
                upstream_exon_tree[es:ee] = 1

        if cds_end_genomic is not None:
            if strand == "+":
                utr_length = max(0, te_end - max(cds_end_genomic, te_start))
            else:
                utr_length = max(0, min(cds_end_genomic, te_end) - te_start)
        else:
            utr_length = te_length

        rec = {
            "gid": gid, "chrom": chrom, "strand": strand,
            "te_start": te_start, "te_end": te_end, "te_length": te_length,
            "cds_end_genomic": cds_end_genomic,
            "utr_length": utr_length,
            "gene_start": gene_start, "gene_end": gene_end,
            "upstream_exon_tree": upstream_exon_tree,
        }
        gene_records[gid] = rec
        if gene_end > gene_start:
            gene_index[chrom][gene_start:gene_end] = gid

    log.info("%s gene records built", f"{len(gene_records):,}")
    return gene_index, gene_records
