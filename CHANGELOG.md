# Changelog

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
