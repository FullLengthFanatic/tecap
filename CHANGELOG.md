# Changelog

## [0.2.0] — unreleased

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
