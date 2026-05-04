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


# Mechanism glossary -- single source of truth for README, HTML report,
# `tecap explain`, and plot captions. Keys are the canonical CATEGORY strings
# from this module; values are dicts with `short`, `what`, and `why`.
MECHANISM_DEFINITIONS = {
    CAPTURED: {
        "short": "Captured",
        "what":  "Read 3' end lands in the terminal exon and the read covers >=50% of it.",
        "why":   "Successful full-length capture of the mRNA 3' end; the goal of any 3'-end protocol.",
    },
    MECH_A_CORRECT: {
        "short": "MechA-correct",
        "what":  "Read 3' end is in the TE 3' UTR within +-25 bp of an annotated PolyASite cluster, but the read covers <50% of the TE.",
        "why":   "Truncated transcript that nonetheless terminates at a real polyA site; common with degraded input or short-fragment library prep.",
    },
    MECH_A_INTERNAL_UTR: {
        "short": "MechA-internalUTR",
        "what":  "Read 3' end is in the TE 3' UTR but not within +-25 bp of any annotated polyA cluster.",
        "why":   "Internal oligo-dT priming on an A-rich stretch inside the UTR; classic mispriming signature.",
    },
    INTERNAL_PRIME_TE_CDS: {
        "short": "IP-TE-CDS",
        "what":  "Read 3' end is inside the terminal exon's CDS portion (i.e. before the UTR begins).",
        "why":   "Internal priming on the coding portion of the TE; strong mispriming signal.",
    },
    MECH_A_NO_CDS: {
        "short": "MechA-noCDS",
        "what":  "Read 3' end is inside the TE of a non-coding gene (no CDS to delimit a UTR).",
        "why":   "Coding/UTR split does not apply; reported separately so the coding-gene buckets stay clean.",
    },
    MECH_B_APA: {
        "short": "MechB-APA",
        "what":  "Read 3' end is upstream of the terminal exon but at an annotated polyA cluster on an upstream exon.",
        "why":   "Alternative polyadenylation isoform; biological, not a mispriming artifact.",
    },
    MECH_B_EXON: {
        "short": "MechB-exon",
        "what":  "Read 3' end is on an upstream exon, with no annotated polyA cluster nearby.",
        "why":   "Internal priming on an upstream exon's A-rich sequence.",
    },
    MECH_B_ASPECI: {
        "short": "MechB-aspecific",
        "what":  "Read 3' end is upstream of the TE in an intron or gene flank.",
        "why":   "Pre-mRNA priming or off-target alignment.",
    },
    MECH_C: {
        "short": "MechC",
        "what":  "Read 3' end is downstream of the terminal exon end.",
        "why":   "Read-through, unannotated 3' UTR extension, or alignment artifact past the annotated TE.",
    },
}


# Same shape as MECHANISM_DEFINITIONS but for the 8 basecomp BUCKETS.
BUCKET_DEFINITIONS = {
    BUCKET_CAPT: {
        "short": "Captured",
        "what":  "Full-length capture: 3' end in the TE, read covers >=50% of the TE.",
        "why":   "Reference for what a real, on-target cleavage downstream window looks like.",
    },
    BUCKET_MA_PAS: {
        "short": "MechA_PAS+",
        "what":  "MechA-correct reads at a polyA cluster carrying a canonical AAUAAA-like hexamer.",
        "why":   "High-confidence real polyadenylation; the downstream window should not be A-rich.",
    },
    BUCKET_MA_NOPAS: {
        "short": "MechA_PAS-",
        "what":  "MechA-correct reads at a polyA cluster without a canonical hexamer.",
        "why":   "Weaker support for real cleavage; downstream-A enrichment here flags mispriming on A-rich UTR.",
    },
    BUCKET_MB_APA_PAS: {
        "short": "MechB_APA_PAS+",
        "what":  "MechB-APA reads at a hexamer-positive upstream polyA cluster.",
        "why":   "Real APA isoform; downstream-A pattern should match Captured.",
    },
    BUCKET_MB_APA_NOPAS: {
        "short": "MechB_APA_PAS-",
        "what":  "MechB-APA reads at an upstream cluster without a hexamer.",
        "why":   "Candidate APA but unsupported; A-tract enrichment here suggests mispriming.",
    },
    BUCKET_MB_EXON: {
        "short": "MechB_exon",
        "what":  "Reads with 3' end on an upstream exon, no polyA cluster.",
        "why":   "Upstream-exon mispriming; the downstream window should be A-rich if oligo-dT-driven.",
    },
    BUCKET_MB_ASPECI: {
        "short": "MechB_aspecific",
        "what":  "Reads with 3' end in an intron or gene flank upstream of the TE.",
        "why":   "Pre-mRNA priming; tends to be A-rich at the cleavage point.",
    },
    BUCKET_IP_CDS: {
        "short": "IP_TE_CDS",
        "what":  "Reads with 3' end in the terminal exon's CDS portion.",
        "why":   "Internal priming in CDS; strong mispriming signal.",
    },
}


# Captions printed under each figure (drawn via fig.text). One sentence,
# dense; the HTML report renders the same text under each PNG.
PLOT_CAPTIONS = {
    "single": (
        "Mechanism breakdown for {sample}. Left panel: fraction of multi-exon "
        "reads in each of the 9 buckets. Centre panel: read-length distribution "
        "for Captured vs MechA-correct. Right panel: rates of Captured and "
        "MechA-correct stratified by 3' UTR length, where mispriming bias "
        "concentrates."
    ),
    "mecha_scatter": (
        "MechA-correct reads only: read length vs fraction of the terminal exon "
        "covered, hexbin density. Reads above the dashed coverage threshold get "
        "promoted to Captured."
    ),
    "comparison_classify": (
        "Cross-sample mechanism comparison, three panels left to right. "
        "Left: every mechanism bucket as a horizontal row, one bar per "
        "sample within each row, % of multi-exon reads. "
        "Middle: per-sample MechA-correct rate (% of reads in the bin "
        "classified MechA-correct), one polyline per sample, x-axis is "
        "the 3' UTR length bin. "
        "Right: per-sample Captured rate by 3' UTR length bin (same axes "
        "as middle, but for Captured)."
    ),
    "basecomp": (
        "Per-bucket histograms of %A in the {window} nt reference window "
        "immediately downstream of each cleavage site. Grey shading marks the "
        "moderate-A regime (30-50% A); the dashed line marks classical A-tract "
        "priming (>=60% A). Enrichment in the grey band reflects priming on "
        "moderate-A sequence; enrichment past the dashed line reflects "
        "priming on classical A-tracts. Empirically, single-cell prep "
        "datasets (10x, BD Rhapsody, ArgenTag, plate FLASH-seq) cluster in "
        "the grey band; bulk Iso-Seq datasets cluster past the dashed line. "
        "The biochemical driver of this split is currently uncharacterized."
    ),
    "comparison_basecomp": (
        "Cross-sample %A histograms in the {window} nt window downstream of "
        "cleavage. Grey band: 30-50% A (moderate-A priming). Dashed line: "
        ">=60% A (classical A-tract). Compare buckets across samples to spot "
        "protocol-specific mispriming."
    ),
}
