"""Per-read classification into the nine mechanism buckets."""

from tecap.constants import (
    CAPTURED, MECH_A_CORRECT, MECH_A_INTERNAL_UTR, INTERNAL_PRIME_TE_CDS,
    MECH_A_NO_CDS, MECH_B_APA, MECH_B_EXON, MECH_B_ASPECI, MECH_C,
)
from tecap.polya import hits_polya, in_upstream_exon


def classify_read(read, gene_rec, polya_index, cov_threshold):
    """Returns (category, cov_frac, read_3prime, has_pas).

    has_pas is True iff the 3' end falls within a polyA cluster that carries
    a canonical PAS hexamer. Only meaningful for MECH_A_CORRECT and MECH_B_APA.
    """
    strand   = gene_rec["strand"]
    te_start = gene_rec["te_start"]
    te_end   = gene_rec["te_end"]
    te_len   = gene_rec["te_length"]
    chrom    = gene_rec["chrom"]
    cds_end  = gene_rec["cds_end_genomic"]

    if strand == "+":
        read_3prime = read.reference_end
    else:
        read_3prime = read.reference_start

    overlap  = max(0, min(read.reference_end, te_end) -
                      max(read.reference_start, te_start))
    cov_frac = overlap / te_len if te_len > 0 else 0.0

    if cov_frac >= cov_threshold:
        return CAPTURED, cov_frac, read_3prime, False

    if strand == "+":
        upstream_of_te   = read_3prime <= te_start
        downstream_of_te = read_3prime >  te_end
    else:
        upstream_of_te   = read_3prime >= te_end
        downstream_of_te = read_3prime <  te_start

    if upstream_of_te:
        hit, has_pas = hits_polya(chrom, read_3prime, strand, polya_index)
        if hit:
            return MECH_B_APA, cov_frac, read_3prime, has_pas
        if in_upstream_exon(read_3prime, gene_rec):
            return MECH_B_EXON, cov_frac, read_3prime, False
        return MECH_B_ASPECI, cov_frac, read_3prime, False

    if downstream_of_te:
        return MECH_C, cov_frac, read_3prime, False

    if cds_end is None:
        return MECH_A_NO_CDS, cov_frac, read_3prime, False

    if strand == "+":
        in_utr = read_3prime > cds_end
    else:
        in_utr = read_3prime < cds_end

    if in_utr:
        hit, has_pas = hits_polya(chrom, read_3prime, strand, polya_index)
        if hit:
            return MECH_A_CORRECT, cov_frac, read_3prime, has_pas
        return MECH_A_INTERNAL_UTR, cov_frac, read_3prime, False
    return INTERNAL_PRIME_TE_CDS, cov_frac, read_3prime, False
