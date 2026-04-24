"""Category, bucket, and plotting constants shared across tecap modules."""

import numpy as np


CAPTURED              = "Captured (>=50%)"
MECH_A_CORRECT        = "MechA-correct: in TE UTR, at polyA site, too short"
MECH_A_INTERNAL_UTR   = "MechA-internalUTR: in TE UTR, no polyA site"
INTERNAL_PRIME_TE_CDS = "Internal priming in TE CDS (no UTR)"
MECH_A_NO_CDS         = "MechA: in TE, non-coding gene"
MECH_B_APA            = "MechB-APA: upstream TE, known polyA site"
MECH_B_EXON           = "MechB-exon: upstream TE, priming on upstream exon"
MECH_B_ASPECI         = "MechB-aspecific: upstream TE, intronic/flanking"
MECH_C                = "MechC: downstream of TE"

CATEGORIES = [
    CAPTURED, MECH_A_CORRECT, MECH_A_INTERNAL_UTR, INTERNAL_PRIME_TE_CDS,
    MECH_A_NO_CDS, MECH_B_APA, MECH_B_EXON, MECH_B_ASPECI, MECH_C,
]

LABELS_SHORT = {
    CAPTURED:              "Captured\n(>=50%)",
    MECH_A_CORRECT:        "MechA-correct:\nin UTR, at polyA",
    MECH_A_INTERNAL_UTR:   "MechA-internalUTR:\nin UTR, no polyA",
    INTERNAL_PRIME_TE_CDS: "Internal priming\nin TE CDS",
    MECH_A_NO_CDS:         "MechA: in TE\nnon-coding",
    MECH_B_APA:            "MechB-APA:\nupstream, known polyA",
    MECH_B_EXON:           "MechB-exon:\nupstream exon priming",
    MECH_B_ASPECI:         "MechB-aspecific:\nintronic/flanking",
    MECH_C:                "MechC:\ndownstream",
}

COLORS = {
    CAPTURED:              "#2ecc71",
    MECH_A_CORRECT:        "#27ae60",
    MECH_A_INTERNAL_UTR:   "#e67e22",
    INTERNAL_PRIME_TE_CDS: "#e74c3c",
    MECH_A_NO_CDS:         "#f39c12",
    MECH_B_APA:            "#3498db",
    MECH_B_EXON:           "#c0392b",
    MECH_B_ASPECI:         "#8e44ad",
    MECH_C:                "#95a5a6",
}

UTR_BINS       = [0, 500, 1000, 2000, 5000, np.inf]
UTR_BIN_LABELS = ["0-500", "500-1k", "1k-2k", "2k-5k", ">5k"]


# Base composition buckets (basecomp subcommand)
BUCKET_CAPT         = "Captured"
BUCKET_MA_PAS       = "MechA_PAS+"
BUCKET_MA_NOPAS     = "MechA_PAS-"
BUCKET_MB_APA_PAS   = "MechB_APA_PAS+"
BUCKET_MB_APA_NOPAS = "MechB_APA_PAS-"
BUCKET_MB_EXON      = "MechB_exon"
BUCKET_MB_ASPECI    = "MechB_aspecific"
BUCKET_IP_CDS       = "IP_TE_CDS"

BUCKETS = [
    BUCKET_CAPT,
    BUCKET_MA_PAS, BUCKET_MA_NOPAS,
    BUCKET_MB_APA_PAS, BUCKET_MB_APA_NOPAS,
    BUCKET_MB_EXON, BUCKET_MB_ASPECI,
    BUCKET_IP_CDS,
]

BUCKET_COLORS = {
    BUCKET_CAPT:         "#2ecc71",
    BUCKET_MA_PAS:       "#27ae60",
    BUCKET_MA_NOPAS:     "#e67e22",
    BUCKET_MB_APA_PAS:   "#3498db",
    BUCKET_MB_APA_NOPAS: "#9b59b6",
    BUCKET_MB_EXON:      "#c0392b",
    BUCKET_MB_ASPECI:    "#8e44ad",
    BUCKET_IP_CDS:       "#e74c3c",
}


SCHEMA_VERSION = "0.1"
