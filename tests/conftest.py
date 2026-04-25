"""Synthetic fixtures: one contig, three genes, a polyA atlas, and a BAM with
reads engineered to land in every classification bucket."""

import os

import pysam
import pytest

CONTIG      = "chr_test"
CONTIG_LEN  = 25000


# GENE_A: + strand, coding. Exons: [1000,1500], [3000,3500], [7000,10000] (TE).
# CDS: [1200..8000]. te_length = 3000, utr_length = 2000.
GENE_A = {
    "gid": "GA",
    "strand": "+",
    "exons": [(1000, 1500), (3000, 3500), (7000, 10000)],
    "cds":   [(1200, 1500), (3000, 3500), (7000, 8000)],
}

# GENE_B: - strand, coding. Exons: [12000,14000] (TE on minus), [14500,15500],
# [17000,18000]. CDS: [12500..17500]. te_start=12000, te_end=14000, te_length=2000.
# cds_end_genomic on - = min CDS start = 12500 → utr_length = 500.
GENE_B = {
    "gid": "GB",
    "strand": "-",
    "exons": [(12000, 14000), (14500, 15500), (17000, 18000)],
    "cds":   [(12500, 14000), (14500, 15500), (17000, 17500)],
}

# GENE_C: + strand, non-coding. Exons: [20000,20500], [22000,24000] (TE).
GENE_C = {
    "gid": "GC",
    "strand": "+",
    "exons": [(20000, 20500), (22000, 24000)],
    "cds":   [],
}

GENES = [GENE_A, GENE_B, GENE_C]


# (chrom, start, end, name, score, strand, col7, col8, col9, type, hexamer)
POLYA_CLUSTERS = [
    (CONTIG,  9900, 10000, "A_te",       "0", "+", "0", "0", "0", "TE", "AATAAA"),
    (CONTIG,  8400,  8500, "A_te_nopas", "0", "+", "0", "0", "0", "TE", "NaN"),
    (CONTIG,  3400,  3500, "A_upstream", "0", "+", "0", "0", "0", "DS", "AATAAA"),
    (CONTIG, 12000, 12100, "B_te",       "0", "-", "0", "0", "0", "TE", "AATAAA"),
    (CONTIG, 17900, 18000, "B_upstream", "0", "-", "0", "0", "0", "DS", "AATAAA"),
]


# (name, gid, start, cigar_ops)  — CIGAR ops: 0=M, 3=N
# reference_start/reference_end forms a box used by classify_read; the CIGAR
# just needs at least one N so multi-exon filter passes.
READS = [
    # GENE_A (TE=[7000,10000], CDS ends at 8000, threshold 0.5):
    # CAPTURED: cov = 2000/3000 = 0.67
    ("r_captured_1",      "GA", 3000, [(0, 500), (3, 3500), (0, 2000)]),
    # MECH_A_CORRECT PAS+: ref=[8600,9950], cov=1350/3000=0.45, 3'=9950
    ("r_mechA_pas_1",     "GA", 8600, [(0, 100), (3, 1200), (0,   50)]),
    # MECH_A_CORRECT PAS-: ref=[7100,8450], cov=1350/3000=0.45, 3'=8450
    ("r_mechA_nopas_1",   "GA", 7100, [(0, 100), (3, 1200), (0,   50)]),
    # MECH_A_INTERNAL_UTR: ref=[7100,8300], cov=1200/3000=0.40, 3'=8300
    ("r_mechA_intutr_1",  "GA", 7100, [(0, 100), (3, 1050), (0,   50)]),
    # INTERNAL_PRIME_TE_CDS: ref=[7100,7500], cov=400/3000=0.13, 3'=7500<8000
    ("r_ip_te_cds_1",     "GA", 7100, [(0, 100), (3,  200), (0,  100)]),
    # MECH_B_APA PAS+: 3'=3500 in [3400,3500] cluster, cov=0
    ("r_mechB_apa_1",     "GA", 1000, [(0, 200), (3, 2100), (0,  200)]),
    # MECH_B_EXON: 3'=3200 in exon2, cov=0
    ("r_mechB_exon_1",    "GA", 1000, [(0, 200), (3, 1800), (0,  200)]),
    # MECH_B_ASPECI: 3'=4200 intronic, cov=0
    ("r_mechB_aspeci_1",  "GA", 1000, [(0, 200), (3, 2800), (0,  200)]),
    # MECH_C: 3'=10500 past te_end
    ("r_mechC_1",         "GA", 9500, [(0, 200), (3,  300), (0,  500)]),

    # GENE_C: non-coding TE, cov=0.1, 3'=22200 in TE → MECH_A_NO_CDS
    ("r_mechA_nocds_1",   "GC", 20200, [(0, 200), (3, 1600), (0, 200)]),

    # GENE_B (minus-strand orientation sanity): CAPTURED. ref=[12000,15100],
    # overlap with TE=[12000,14000] = 2000/2000 = 1.0.
    ("r_captured_B",      "GB", 12000, [(0, 200), (3, 2800), (0, 100)]),
]


# The bucket each read is expected to land in (used by test_classify + test_bam).
EXPECTED_BUCKETS = {
    "r_captured_1":     "Captured (>=50%)",
    "r_mechA_pas_1":    "MechA-correct: in TE UTR, at polyA site, too short",
    "r_mechA_nopas_1":  "MechA-correct: in TE UTR, at polyA site, too short",
    "r_mechA_intutr_1": "MechA-internalUTR: in TE UTR, no polyA site",
    "r_ip_te_cds_1":    "Internal priming in TE CDS (no UTR)",
    "r_mechB_apa_1":    "MechB-APA: upstream TE, known polyA site",
    "r_mechB_exon_1":   "MechB-exon: upstream TE, priming on upstream exon",
    "r_mechB_aspeci_1": "MechB-aspecific: upstream TE, intronic/flanking",
    "r_mechC_1":        "MechC: downstream of TE",
    "r_mechA_nocds_1":  "MechA: in TE, non-coding gene",
    "r_captured_B":     "Captured (>=50%)",
}

# PAS+/- expectation for MECH_A_CORRECT / MECH_B_APA reads.
EXPECTED_HAS_PAS = {
    "r_mechA_pas_1":   True,
    "r_mechA_nopas_1": False,
    "r_mechB_apa_1":   True,
}


def _build_fasta(path):
    seq = bytearray(b"C" * CONTIG_LEN)
    # Planted A-stretches for basecomp tests:
    seq[9950:9970] = b"A" * 20   # downstream of r_mechA_pas_1's 3' end
    seq[8450:8470] = b"A" * 20   # downstream of r_mechA_nopas_1's 3' end

    with open(path, "w") as f:
        f.write(f">{CONTIG}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60].decode() + "\n")
    pysam.faidx(path)
    return path


def _build_gtf(path):
    lines = []
    for g in GENES:
        attrs = f'gene_id "{g["gid"]}"; gene_name "{g["gid"]}_sym";'
        for (s, e) in g["exons"]:
            lines.append(f"{CONTIG}\thavana\texon\t{s+1}\t{e}\t.\t{g['strand']}\t.\t{attrs}")
        for (s, e) in g["cds"]:
            lines.append(f"{CONTIG}\thavana\tCDS\t{s+1}\t{e}\t.\t{g['strand']}\t0\t{attrs}")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return path


def _build_polya_bed(path):
    with open(path, "w") as f:
        for row in POLYA_CLUSTERS:
            f.write("\t".join(str(x) for x in row) + "\n")
    return path


def _build_bam(sorted_path):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}],
    }
    unsorted = sorted_path + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as bam:
        for (name, gid, start, cigar) in READS:
            a = pysam.AlignedSegment()
            a.query_name = name
            qlen = sum(length for op, length in cigar if op == 0)
            a.query_sequence  = "A" * qlen
            a.query_qualities = pysam.qualitystring_to_array("I" * qlen)
            a.flag            = 16 if gid == "GB" else 0  # reverse for minus gene
            a.reference_id    = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigartuples     = cigar
            if gid == "GA":
                a.set_tag("CB", f"CELL_{name[-1]}")
            bam.write(a)
    pysam.sort("-o", sorted_path, unsorted)
    pysam.index(sorted_path)
    os.remove(unsorted)
    return sorted_path


@pytest.fixture(scope="session")
def synthetic_fixture(tmp_path_factory):
    d = tmp_path_factory.mktemp("tecap_fixture")
    gtf_path   = _build_gtf(str(d / "test.gtf"))
    polya_path = _build_polya_bed(str(d / "polya.bed"))
    fasta_path = _build_fasta(str(d / "chr_test.fa"))
    bam_path   = _build_bam(str(d / "test.bam"))
    return {
        "dir":   str(d),
        "gtf":   gtf_path,
        "polya": polya_path,
        "fasta": fasta_path,
        "bam":   bam_path,
        "contig": CONTIG,
        "contig_len": CONTIG_LEN,
        "reads":    READS,
        "expected_buckets": EXPECTED_BUCKETS,
        "expected_has_pas": EXPECTED_HAS_PAS,
    }
