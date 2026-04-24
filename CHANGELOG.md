# Changelog

## [0.1.0] — unreleased

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
