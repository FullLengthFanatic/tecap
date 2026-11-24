# Changelog

## [0.5.0] - <date of release, set by Simone>

tecap now puts a read's 3' end where flight-polyafilter-bench 1.2.0 puts
its priming site, and labels each site templated, non-templated or
ambiguous. The default site moves from the alignment end to the start of
the read's terminal adenine run, and on stranded input a read antisense
to its best-overlap gene is classified against a gene on its own strand,
or skipped when there is none, so most numbers of `classify` and `basecomp`
change (old and new values under Changed). `--site-anchor alignment_end
--strandedness unstranded` gives the 0.4.0 numbers back, apart from the
off-by-one under Fixed. JSON schema 0.2.

### Added
- `--site-anchor {alignment_end,a_run_start}` on `classify` and
  `basecomp`, default `a_run_start`. The anchor sets the read 3' end used
  in every position test of `classify` and, in `basecomp`, where the
  downstream window starts. `a_run_start` is the site rule of
  flight-polyafilter-bench 1.2.0, ported in the new module `tecap.site`:
  at the gene-strand 3' end, the read's terminal adenine run is grown over
  the soft clip and the aligned bases with one interior non-A allowed, and
  when it holds 10 or more adenines the site is the reference position
  where the run starts (the alignment end when the run lies wholly in the
  clip). A read without such a run, and a record stored without SEQ, falls
  back to the alignment end. Coverage of the terminal exon stays on the
  alignment, so the Captured call does not depend on the anchor.
- `site_anchor` and `site_counts` in the `classify` and `basecomp` JSON.
  `site_counts` is `{"observed": n, "fallback": n}`: reads whose site came
  from an adenine run, and reads that fell back to the alignment end.
  Under `alignment_end` every read counts as observed.
- Labels in the `basecomp` JSON. A read with an observed site is labelled
  from the integer count of adenines in the 20 nt of reference downstream
  of its site: templated (8 or more), non-templated (4 or fewer) or
  ambiguous (5 to 7). New keys per bucket: `n_site`, `n_templated`,
  `n_nontemplated`, `n_ambiguous`, `frac_templated`, `frac_nontemplated`,
  `frac_ambiguous`. The fractions are over `n_site`, not `n`, and are null
  when `n_site` is 0. With `--window` other than 20 the seven keys are null
  and a warning is logged; the bands are computed as before.
- New `thresholds` keys in the `basecomp` JSON: `templated_min` (8),
  `nontemplated_max` (4), `label_window` (20), `site_run_min` (10),
  `site_run_tolerance` (1).
- `definition_source` in the `basecomp` JSON
  (`"flight-polyafilter-bench 1.2.0"`), the definition the site rule and
  the labels come from.
- Bucket by label table in the HTML report (`tecap report`, single sample
  and comparison): per bucket, the reads with a site and how many of them
  are templated, non-templated and ambiguous, with a caption. Left out for
  a basecomp JSON without labels.
- Stranding guard. `classify` and `basecomp` now emit a loud warning when
  the BAM/gene-strand match fraction falls below `ORIENT_MIN_FRAC` (0.80,
  in `constants.py`), which flags unstranded input (e.g. ONT cDNA not run
  through pychopper / cdna_classifier). The read 3' end and the basecomp
  downstream window are defined from the gene strand, so unstranded
  input silently dilutes those metrics toward the genomic mean. `basecomp`
  now computes orientation itself and records `orient_match_frac` in its
  JSON (previously only `classify` did).
- `--strandedness {auto,stranded,unstranded}` on `classify` and
  `basecomp`, default `auto`. A read's gene is the gene whose span
  overlaps it the most, genes on both strands competing. On stranded
  input a read antisense to that gene (read strand opposite the gene
  strand) is classified against the best-overlap gene on its own strand
  instead, and gets the category, site, window, per-cell and per-gene
  count that gene gives; a read with no gene on its strand is skipped:
  it gets no category, site, window, per-cell or per-gene count. On
  unstranded input every read keeps its best-overlap gene on either
  strand. `auto` reads the first 200,000 reads the run would classify,
  contigs in BAM order, and calls the input stranded when 0.80
  (`ORIENT_MIN_FRAC`) or more of them match the strand of their
  best-overlap gene on either strand. Three new keys in both JSON
  outputs: `strandedness`, `{"requested", "resolved", "probe"}`, where
  `probe` is `{"reads", "match_frac", "stranded"}` under `auto` and null
  otherwise; `antisense_reassigned`, the number of reads classified
  against a gene on their own strand; and `antisense_skipped`, the number
  of reads skipped. `orient_match_frac` is computed against the
  best-overlap gene on either strand, before the reassignment and the
  skip, and keeps its meaning. Python API: `bam.gene_of_read` (its
  `strand` argument lets only genes of one strand compete),
  `bam.choose_gene`, which both workers use for the gene and the
  antisense status, `bam.probe_strandedness`, and `skip_antisense` on
  `bam.analyse_bam` and `basecomp.analyse`, which count the reassigned
  and skipped reads under `antisense_reassigned` and `antisense_skipped`
  (in the `analyse_bam` result and in the `site_stats` of `analyse`).
- Base composition thresholds stamped into the `basecomp` JSON. The cutoffs
  that produce `frac_ge_60` / `frac_30_50` / `frac_ge_40` and the histogram
  bin count are recorded under a top-level `thresholds` key and sourced from
  new `BASECOMP_*` constants, so a basecomp JSON is self-describing rather
  than depending on the code version that wrote it. Additive change, and
  existing consumers are unaffected; the schema moves to 0.2 with the
  other changes of this release (see Changed).
- `whitepaper.md`: a self-contained explainer of what tecap does, its
  strengths, and its limitations, for readers new to the tool.

### Changed
- Default read 3' end: the site under `a_run_start` (see Added), in place
  of the alignment end. For a read with aligned adenines at its 3' end the
  site lies upstream of the alignment end, so the read can change bucket,
  and its `basecomp` window now starts where its adenine run starts and
  covers the reference adenines the read aligned to. Numbers the anchor
  changes:
  - `basecomp`, per bucket: `n` (every bucket except `Captured`),
    `median`, `mean`, `hist`, `frac_ge_60`, `frac_30_50`, `frac_ge_40`.
  - `classify`: `counts` and `fractions` of every category except
    `Captured (>=50%)`, and the numbers computed per category:
    `pas_counts`, `pas_fraction`, `utr_bin_mecha_correct`,
    `utr_bin_mecha_rate`, `cov_median_by_cat`, `read_len_median_mecha`,
    `per_cell_metrics`, the per-gene TSV, and the MultiQC columns
    `mecha_correct_pct`, `mechb_aspeci_pct`, `mecha_pas_pos_frac` and
    `mechb_apa_pas_pos_frac`.

  `--site-anchor alignment_end --strandedness unstranded` reproduces
  0.4.0, apart from the off-by-one under Fixed. Old and new values are in
  the tables at the end of this section.
- On stranded input a read antisense to its best-overlap gene is
  classified against the best-overlap gene on its own strand, or skipped
  when no gene on its strand overlaps it (`--strandedness`, see Added).
  The default `auto` calls both libraries
  below stranded (probe match fraction 0.984 for BD67 and 0.978 for
  10x_FL on the FLNC BAM, 0.994 and 0.988 on chr21 of `accepted.bam`), so
  their 0.5.0 numbers leave the skipped reads out. Reads skipped, of the reads
  assigned to a gene: `basecomp` on the whole FLNC BAM, 1,755,777 of
  50,538,102 (3.47%) for BD67 and 2,214,195 of 52,814,167 (4.19%) for
  10x_FL; `classify` on chr21 of `accepted.bam`, 6,073 of 438,044 (1.39%)
  and 7,405 of 410,831 (1.80%). In `basecomp` MechB_aspecific loses the
  most: under the default anchor 1,359,294 of 3,323,636 reads (40.9%) for
  BD67 and 1,646,175 of 5,440,400 (30.3%) for 10x_FL, and its bands rise
  (BD67 `frac_ge_60` 0.280 without the skip, 0.471 with it). The other
  buckets lose 0.01% to 1.8% of their reads. Numbers that change on
  stranded input, under either anchor:
  - `basecomp`: `n` of every bucket, `Captured` included, with `hist`,
    `mean` and the three bands; `median` in some buckets; the label keys
    (in every bucket under `alignment_end`, in some under the default
    anchor, as most antisense reads have no observed site); `site_counts`.
  - `classify`: `total`, `site_counts`, and `counts` and `fractions` of
    every category, `Captured (>=50%)` included; `pas_counts`,
    `pas_fraction`, `utr_bin_total`, `utr_bin_captured`,
    `utr_bin_capt_rate`, `utr_bin_mecha_correct`, `utr_bin_mecha_rate`,
    `cov_median_by_cat`, `read_len_median_capt`, `per_cell_metrics`, the
    per-gene TSV (genes with only antisense reads drop out: 247 to 237
    genes for BD67 and 289 to 269 for 10x_FL on chr21), and the MultiQC
    columns `captured_pct`, `mecha_correct_pct`, `mechb_aspeci_pct`,
    `mecha_pas_pos_frac`, `mechb_apa_pas_pos_frac` and
    `total_classified`. `read_len_median_mecha` did not change on chr21
    but can.

  `orient_match_frac` does not change: it is computed before the
  reassignment and the skip. With `--strandedness unstranded`, or under
  `auto` on input it calls unstranded, nothing is reassigned or skipped.
- `schema_version` is 0.2 in both JSON outputs: new keys (see Added), and
  under the default anchor the numbers are not comparable with those of a
  schema 0.1 JSON.
- Wording. The adenine bands are descriptive and no longer name a priming
  mechanism: "moderate-A priming" and "classical A-tract priming" are gone
  from the README, the whitepaper, the HTML report, the plot captions, the
  `basecomp` plot legend (now "30% to 50% A" and "60% A or more"),
  `CITATION.cff` and `.zenodo.json`. The bucket explanations
  (`tecap explain`, the report legend) and the MultiQC MechB-aspecific
  column description no longer infer a priming mechanism from a bucket. A
  bucket says where the site is; the label says what kind of site it is.
  The report glossary gains `site` and `label` rows.
- The entries for 0.4.0 and earlier below keep the retired wording by
  design.
- The example plots and numbers in `README.md` and `whitepaper.md` come
  from runs of 0.4.0 or earlier, at the alignment end, and were not
  regenerated for this release.
- Python API: `basecomp.analyse` returns
  `(a_count_by_bucket, bucket_counts, orient, site_stats)`, integer
  adenine counts in place of shares plus a fourth element;
  `basecomp.summarise` takes those counts and `window` and `labels`;
  `bam.analyse_bam` returns a `site` key; `classify.classify_read` takes
  `read_3prime`. `basecomp.pct_a` is removed.
- Disambiguated `BD46_FS_SEQ` in the README and basecomp plot
  caption: it is an in-house FLASH-seq variant adapted to the BD
  Rhapsody cartridge/beads with Kinnex concatenation, not stock
  BD Rhapsody chemistry. Earlier wording ("BD Rhapsody Kinnex")
  conflated a custom protocol with the off-the-shelf platform.

Old and new values for two in-house PacBio HiFi libraries, `BD67`
(FLIGHT-seq) and `10x_FL` (10x); they are not a matched pair. References:
Cell Ranger GRCh38-2024-A FASTA and GTF, PolyASite 3.0 atlas
(`atlas.clusters.3.0.GRCh38.GENCODE_42`), `--polya-slop 25`,
`--platform cdna-pacbio`. The 0.5.0 runs use the default `--strandedness
auto`, which calls both libraries stranded, so their antisense reads are
skipped.

`basecomp` on the whole FLNC BAM, `--window 20`. `n` is given for 0.4.0,
for 0.5.0 with `--site-anchor alignment_end` (the Fixed entry and the
antisense skip) and for 0.5.0 with the default; the other cells read
0.4.0 / 0.5.0 default.

BD67:

| bucket | `n` 0.4.0 | `n` `alignment_end` | `n` 0.5.0 | `median` | `mean` | `frac_ge_60` | `frac_30_50` | `frac_ge_40` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Captured | 22,204,457 | 22,042,556 | 22,042,556 | 0.15 / 0.20 | 0.154 / 0.250 | 0.005 / 0.076 | 0.118 / 0.213 | 0.047 / 0.182 |
| MechA_PAS+ | 7,369,742 | 7,360,616 | 7,344,817 | 0.20 / 0.25 | 0.201 / 0.312 | 0.010 / 0.128 | 0.208 / 0.304 | 0.097 / 0.283 |
| MechA_PAS- | 2,143,146 | 2,125,640 | 2,107,229 | 0.30 / 0.70 | 0.278 / 0.616 | 0.023 / 0.619 | 0.413 / 0.143 | 0.241 / 0.792 |
| MechB_APA_PAS+ | 2,058,349 | 2,035,295 | 1,953,001 | 0.40 / 0.55 | 0.348 / 0.470 | 0.101 / 0.453 | 0.391 / 0.152 | 0.517 / 0.692 |
| MechB_APA_PAS- | 3,809,570 | 3,788,237 | 3,974,590 | 0.35 / 0.50 | 0.366 / 0.494 | 0.048 / 0.372 | 0.527 / 0.306 | 0.475 / 0.740 |
| MechB_exon | 7,280,915 | 7,595,057 | 7,583,909 | 0.25 / 0.30 | 0.246 / 0.377 | 0.021 / 0.220 | 0.302 / 0.273 | 0.215 / 0.421 |
| MechB_aspecific | 3,807,115 | 2,052,650 | 1,964,342 | 0.20 / 0.55 | 0.233 / 0.551 | 0.019 / 0.471 | 0.302 / 0.236 | 0.185 / 0.645 |
| IP_TE_CDS | 1,288,723 | 1,283,380 | 1,314,216 | 0.35 / 0.55 | 0.352 / 0.521 | 0.064 / 0.446 | 0.429 / 0.247 | 0.471 / 0.760 |

10x_FL:

| bucket | `n` 0.4.0 | `n` `alignment_end` | `n` 0.5.0 | `median` | `mean` | `frac_ge_60` | `frac_30_50` | `frac_ge_40` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Captured | 17,050,085 | 16,854,141 | 16,854,141 | 0.15 / 0.20 | 0.168 / 0.248 | 0.005 / 0.078 | 0.146 / 0.211 | 0.064 / 0.190 |
| MechA_PAS+ | 6,903,399 | 6,891,416 | 6,819,914 | 0.20 / 0.25 | 0.206 / 0.306 | 0.008 / 0.131 | 0.220 / 0.267 | 0.110 / 0.281 |
| MechA_PAS- | 2,156,472 | 2,126,324 | 2,139,418 | 0.30 / 0.60 | 0.278 / 0.589 | 0.021 / 0.552 | 0.419 / 0.180 | 0.241 / 0.829 |
| MechB_APA_PAS+ | 2,848,163 | 2,803,405 | 2,724,707 | 0.40 / 0.55 | 0.371 / 0.509 | 0.099 / 0.485 | 0.445 / 0.150 | 0.549 / 0.791 |
| MechB_APA_PAS- | 6,592,149 | 6,551,274 | 6,804,851 | 0.35 / 0.55 | 0.362 / 0.511 | 0.048 / 0.373 | 0.515 / 0.273 | 0.463 / 0.826 |
| MechB_exon | 9,048,790 | 9,187,025 | 9,105,821 | 0.25 / 0.35 | 0.264 / 0.393 | 0.022 / 0.235 | 0.337 / 0.270 | 0.246 / 0.497 |
| MechB_aspecific | 5,801,012 | 3,883,592 | 3,794,225 | 0.25 / 0.60 | 0.275 / 0.590 | 0.031 / 0.568 | 0.367 / 0.178 | 0.264 / 0.807 |
| IP_TE_CDS | 1,768,090 | 1,753,897 | 1,811,570 | 0.35 / 0.55 | 0.348 / 0.528 | 0.057 / 0.428 | 0.420 / 0.260 | 0.461 / 0.835 |

Under the default, `site_counts` is 43,780,242 observed and 5,002,083
fallback for BD67, 41,339,290 and 9,260,682 for 10x_FL.

`classify` on chr21 of the reads that flight-polyafilter-bench stage 30
accepts and passes to tecap (`<lib>/tecap/accepted.bam`), with its options
`--cell-barcode-tag XK --per-gene-table`. Category counts:

| category | BD67 0.4.0 | BD67 `alignment_end` | BD67 0.5.0 | 10x_FL 0.4.0 | 10x_FL `alignment_end` | 10x_FL 0.5.0 |
|---|---:|---:|---:|---:|---:|---:|
| Captured (>=50%) | 200,666 | 196,518 | 196,518 | 105,086 | 101,810 | 101,810 |
| MechA-correct | 56,326 | 56,319 | 56,317 | 48,043 | 47,997 | 47,993 |
| MechA-internalUTR | 467 | 467 | 467 | 451 | 439 | 439 |
| Internal priming in TE CDS (no UTR) | 9,880 | 9,877 | 9,879 | 19,136 | 19,122 | 19,126 |
| MechA: in TE, non-coding gene | 1,124 | 1,067 | 1,067 | 1,116 | 1,055 | 1,055 |
| MechB-APA | 57,993 | 57,921 | 58,665 | 104,811 | 104,587 | 105,906 |
| MechB-exon | 91,693 | 92,829 | 93,228 | 90,482 | 90,662 | 90,508 |
| MechB-aspecific | 19,885 | 16,968 | 15,825 | 41,648 | 37,708 | 36,543 |
| MechC | 10 | 5 | 5 | 58 | 46 | 46 |
| total | 438,044 | 431,971 | 431,971 | 410,831 | 403,426 | 403,426 |

Between `alignment_end` and the default, 1,960 of 431,971 BD67 reads
(0.45%) and 2,676 of 403,426 10x_FL reads (0.66%) change category, most
of them from MechB-aspecific to MechB-exon (1,161 and 1,233) and from
MechB-exon to MechB-APA (757 and 1,332). MechB-APA `pas_fraction` is
0.372 in 0.4.0, 0.371 under `alignment_end` and 0.259 under the default
for BD67; 0.322, 0.313 and 0.263 for 10x_FL.

### Fixed
- Plus-strand off-by-one in `classify`, which also assigns the `basecomp`
  buckets. On a plus-strand gene the read 3' end is the first reference
  base after the read, and 0.4.0 tested that base, not the last cDNA
  base, against the upstream exons and the polyA cluster windows. A read
  ending on the last base of an upstream exon was called MechB-aspecific
  instead of MechB-exon, and a read ending one base before an exon or a
  cluster window was called MechB-exon or MechB-APA instead of
  MechB-aspecific. The 3' UTR polyA test (MechA-correct against
  MechA-internalUTR) and the PAS+ / PAS- split had the same shift.
  Minus-strand genes were right. The tests now use the last cDNA base
  under both anchors. What it moves, 0.4.0 against 0.5.0 with
  `--site-anchor alignment_end --strandedness unstranded`, so that the
  antisense skip does not enter:
  - `classify`, chr21: MechB-aspecific 19,885 to 18,348 and MechB-exon
    91,693 to 93,230 (BD67); MechB-aspecific 41,648 to 40,937, MechB-exon
    90,482 to 91,194 and MechB-APA 104,811 to 104,810 (10x_FL).
    `pas_fraction` of MechB-APA 0.372 to 0.371 (BD67) and 0.322 to 0.313
    (10x_FL), of MechA-correct 0.844 to 0.843 and 0.800 to 0.799.
  - `basecomp` `n`: MechB_exon +396,104 and MechB_aspecific -395,171
    (BD67), +276,590 and -271,242 (10x_FL), and smaller changes in the
    PAS+ and PAS- buckets of MechA and MechB-APA (7,514 to 40,803 reads
    per bucket). `median` of BD67 MechB_exon 0.25 to 0.20 and of BD67
    MechB_aspecific 0.20 to 0.25. In the six buckets whose `n` moves,
    `hist` changes and `mean` and the three bands shift by 0.023 at most
    (BD67 MechB_aspecific `frac_30_50`). Captured and IP_TE_CDS do not
    move; their `mean` differs in the seventh significant digit, from
    float32 summation over contigs merged in completion order, which can
    differ between any two runs with `--threads` above 1.

## [0.4.0] — 2026-05-07

### Added
- Self-explanatory prose blocks in HTML reports for non-expert
  readers (manuscript reviewers, conference attendees, collaborators
  who haven't run tecap themselves). Three additions, sourced from
  new constants in `tecap/constants.py`:
  - `REPORT_INTRO`: top-level intro paragraph rendered above every
    single-sample and multi-sample report. Names what the report is,
    what was measured, and where the input came from.
  - `HOW_TO_READ_COMPARE`: reading guide rendered only in
    multi-sample comparison reports. Maps polylines / grouped bars /
    colours back to samples.
  - `GLOSSARY`: two-column glossary table at the bottom of every
    report. Defines TE, UTR, CDS, polyA, PAS, APA, oligo-dT.
- New report helpers `_intro_html`, `_how_to_read_html`,
  `_glossary_html` in `tecap/report.py`. Wired into both
  `build_single_report` and `build_compare_report`.
- New CSS classes `.intro` and `.howto` for visually framed prose
  callouts (matches existing `.tile` styling).

## [0.3.2] — 2026-05-04

Prose-only patch. No code, schema, CLI, or behavior change.

### Changed
- Removed an unsupported mechanistic attribution from the README,
  CITATION.cff, .zenodo.json, and the `basecomp` plot caption emitted
  by `tecap/constants.py`. Earlier wording attributed the moderate-A
  vs classical-A split to "saturating-local-concentration oligo-dT
  chemistries (10x GEM droplets, BD Rhapsody capture beads)". That
  claim is not supported by available bench data: bulk Iso-Seq final
  oligo-dT is 1.2 µM, in line with FS-ONT and BD Rhapsody, so bulk
  oligo-dT concentration cannot be the axis driving the split.
- All retracted prose has been replaced with empirical-only language:
  single-cell prep datasets (10x, BD Rhapsody, ArgenTag, plate
  FLASH-seq) cluster in the 30-50% A regime; bulk Iso-Seq datasets
  cluster past the >=60% A line. The biochemical driver of the split
  is currently uncharacterized.

### Note on [0.3.0] history
The `[0.3.0]` block below still contains the original
saturating-local-concentration prose. That entry is preserved verbatim
because v0.3.0 / v0.3.1 are already published on GitHub, Zenodo, and
bioconda; rewriting locked release notes would be misleading. This
`[0.3.2]` entry documents the retraction.

## [0.3.1] — 2026-05-01

CLI ergonomics. JSON schema unchanged.

### Changed
- `tecap report` now accepts space-separated paths for
  `--classify-json` and `--basecomp-json` (`nargs="+"`) instead of
  comma-separated. Commas are legal in POSIX paths, so the comma form
  was a real ambiguity. The new form matches standard argparse
  conventions.

### Migration
- Replace any comma-separated invocation
  `tecap report --classify-json A.json,B.json --basecomp-json A_bc.json,B_bc.json ...`
  with space-separated
  `tecap report --classify-json A.json B.json --basecomp-json A_bc.json B_bc.json ...`.

## [0.3.0] — 2026-04-30

Readability, reporting, and ergonomics. JSON schema unchanged.

### Added
- Single source of truth for mechanism / bucket prose in
  `tecap/constants.py` (`MECHANISM_DEFINITIONS`, `BUCKET_DEFINITIONS`,
  `PLOT_CAPTIONS`). README, HTML report, plot captions, and `tecap explain`
  all render from these dicts.
- `tecap explain [--mechanism NAME] [--scope classify|basecomp|all]
  [--format text|json]` — print the glossary at the terminal.
- `tecap report --classify-json A.json[,B.json,...]
  [--basecomp-json ...] --out-html OUT.html` — single self-contained HTML
  per sample, or cross-sample comparison. Embeds PNGs as base64. No JS,
  no external CSS, no CDN.
- `tecap classify` / `basecomp` accept `--genome {GRCh38,GRCm38,GRCm39}`
  to auto-fetch missing references via the existing
  `fetch_polya_atlas` / `fetch_gencode_gtf`. Adds `--gtf-version` and
  `--ref-cache` (default `$XDG_CACHE_HOME/tecap` / `~/.cache/tecap`).
  Explicit `--polya-sites` / `--gtf` still win.

### Changed
- `comparison_terminal_exon.png` panel 1 is now horizontal grouped bars,
  one row per category, samples grouped within each row. Long category
  names are readable; the previous overlapped vertical x-axis labels are
  gone.
- `comparison_terminal_exon.png` panels 2 and 3 (UTR-bin rates) switched
  from grouped vertical bars to line+marker plots, one polyline per
  sample. Bars hid samples with low MechA-correct rates at N>=4; lines
  scale to any sample count. Caption rewritten to describe left/middle/
  right panels (was "top/bottom").
- `comparison_basecomp.png` switched from side-by-side bars to step-line
  overlay, one polyline per sample per bucket, with a single
  figure-level sample legend. Side-by-side bars were unreadable at
  N>=4. Figure size scales with sample count.
- All four plotting functions now draw a figure-level caption explaining
  what the plot shows.
- Basecomp PNGs now carry an explicit figure-level legend for the grey
  band ("30-50% A: moderate-A priming") and the dashed line
  (">=60% A: classical A-tract"). Per-bucket subplot titles include the
  bucket's interpretation.
- Moderate-A priming attribution updated based on a 4-sample comparison
  (10x Kinnex, BD Rhapsody Kinnex, PacBio Kinnex bulk cerebellum, PacBio
  Kinnex bulk heart, all human GRCh38). MechB_aspecific frac[30,50]:
  10x 0.36, BD46 0.25, Kinnex cerebellum 0.16, Kinnex heart 0.16. Bulk
  Iso-Seq samples instead show heavy classical-A enrichment (frac>=60
  ~ 0.47). The moderate-A signature is characteristic of saturating-
  local-concentration oligo-dT chemistries: 10x GEM droplets (gel bead
  dissolves and releases oligo-dT into a ~1 nL droplet) and BD Rhapsody
  capture beads (oligo-dT density at the bead surface). Free oligo-dT
  at standard concentrations (Iso-Seq, ~20 µL RT) shows classical-A
  internal priming instead. The previous attribution to "saturating
  in-solution oligo-dT" generally was too broad: 10x and Iso-Seq both
  use in-solution oligo-dT, but only 10x's droplet volume creates the
  saturating local concentration that drives moderate-A priming.
  Captions, README, CITATION.cff, and .zenodo.json updated.

### Fixed
- `download-atlas`: PolyASite URLs in `download.py` 404'd against the live
  service. Replaced with the actual paths under
  `polyasite.unibas.ch/download/atlas/{2.0,3.0}/...`. Human is v3.0
  (GENCODE_42), mouse is v2.0 (GRCm38.96) — PolyASite v3.0 mouse is not
  published.
- `--genome` now accepts `{GRCh38, GRCm38, GRCm39}` to reflect the asymmetric
  PolyASite/GENCODE coverage (GRCh38 has both; GRCm38 has PolyASite only;
  GRCm39 has GENCODE only).

### Known limitations / open questions
- The 4-sample chemistry comparison covers droplet-scale (10x), bead-
  surface (BD Rhapsody), and bulk-tube (Kinnex Iso-Seq, ~20 µL RT)
  oligo-dT environments. It does **not** cover plate-scale RT (1-10 µL
  per well, in-solution oligo-dT, e.g. Smart-seq2 / Smart-seq3 /
  FLASH-seq), which would test whether the moderate-A signature tracks
  reaction-volume scale or chemistry lineage. FLASH-seq amplification
  (used in BD Rhapsody) is a Smart-seq descendant, so plate Smart-seq
  could plausibly look like BD46 (chemistry-driven moderate-A) or like
  Iso-Seq (scale-driven, no moderate-A).
- This gap exists because no clean public **human** Smart-seq2/3 +
  PacBio HiFi dataset was findable as of 2026-04-29: PacBio's
  `Kinnex-single-cell-RNA` and `MAS-Seq` buckets are 10x-only;
  HIT-scISOseq corneal limbus is Smart-seq2-derived but not directly
  downloadable; Al'Khafaji TIL T cell data (dbGaP phs003200) is
  10x-derived. Mouse Smart-seq2 + Iso-Seq exists (SRP225196) but
  cross-species adds noise on top of PolyASite v2.0 cluster-type
  filter mismatches.
- Follow-up paths for v0.4: (a) request HIT-scISOseq raw BAMs from
  Zheng/Chen et al. (Sun Yat-sen University); (b) watch PacBio's
  public bucket for plate Smart-seq deposits; (c) generate an in-house
  Smart-seq3 + Kinnex library if the question stays open.

## [0.2.0] — 2026-04-25

### Performance
- Multi-threaded BAM analysis no longer pickles `gene_index`, `gene_records`,
  and `polya_index` per task. The pool is forked from a parent that has
  placed these in module-level globals, so workers inherit them via
  copy-on-write. On a 10x Kinnex BAM, this removes the master-process CPU
  saturation that previously made `--threads N>1` slower than `--threads 1`.
- Contigs with `mapped == 0` in the BAM index are skipped at the parent.
  Avoids ~165 wasted worker dispatches for typical 10x Kinnex BAMs (195
  declared, ~30 populated).

### Added
- `tecap download-atlas --genome {GRCh38,GRCm39}` fetches PolyASite 3.0,
  computes SHA256, and (with `--gtf-version`) optionally fetches a matched
  GENCODE GTF release.
- `{sample}_tecap_mqc.json` MultiQC custom-content table written alongside
  the main classify JSON. Captured/MechA-correct/MechB-aspecific %, PAS+
  fractions, and orientation-match fraction are surfaced.
- `conda/meta.yaml` bioconda recipe stub.
- `docker/Dockerfile` minimal micromamba image.
- `.github/workflows/test.yml` runs pytest on Python 3.10/3.11/3.12 and
  ruff lint on push and PR.
- `CITATION.cff` and `.zenodo.json` so the first Zenodo release picks up
  metadata automatically.

### Unchanged
- JSON output of `classify` and `basecomp` is byte-for-byte identical to
  v0.1 on the same input. No schema bump.

## [0.1.0] — 2026-04-24

First packaged release. Refactor of the two standalone scripts
(`terminal_exon_coverage.py`, `base_composition_downstream.py`) into an
importable package with a proper CLI.

### Added
- `tecap` CLI with subcommands: `classify`, `basecomp`, `compare`.
- `--threads N` multiprocessing by chromosome.
- `--cell-barcode-tag` opt-in per-cell breakdown. Warns and no-ops when the
  tag is absent on the first 10 000 reads.
- `--platform {cdna-pacbio, cdna-ont}` flag. Direct-RNA is refused with a
  clear error message at CLI parse time.
- `--per-gene-table` optional TSV output with per-gene bucket counts.
- JSON schema versioning (`schema_version: "0.1"`).
- Unit tests on a synthetic BAM fixture covering all nine buckets.

### Changed
- `print(...)` replaced with `logging` module. `--verbose` sets DEBUG level.

### Unchanged
- Output JSON schema is backward compatible. Existing consumers work.
- Output filenames unchanged: `{sample}_terminal_exon.{json,png}`, etc.
- All classification logic ported verbatim from the standalone scripts.
