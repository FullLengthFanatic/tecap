# tecap

[![tests](https://github.com/FullLengthFanatic/tecap/actions/workflows/test.yml/badge.svg)](https://github.com/FullLengthFanatic/tecap/actions/workflows/test.yml)

3' terminal exon capture diagnostics for long-read single-cell RNA-seq.

`tecap` classifies long-read alignments by where their 3' end lands relative to the terminal exon (TE), its UTR, and a polyA site atlas. It decomposes capture failures into nine mechanism buckets (successful capture, truncation at a real polyA site, internal priming in the UTR, internal priming in the CDS, alternative polyadenylation, upstream-exon mispriming, intronic mispriming, downstream readthrough) and measures reference base composition downstream of each cleavage site to distinguish classic A-tract internal priming from moderate-A priming characteristic of saturating in-solution oligo-dT.

Designed for PacBio Iso-Seq / Kinnex and Oxford Nanopore cDNA BAMs. Direct-RNA sequencing is explicitly unsupported (no RT, no priming artifact to diagnose).

## Install

```bash
pip install git+https://github.com/FullLengthFanatic/tecap@v0.1.0
```

Development install:

```bash
git clone https://github.com/FullLengthFanatic/tecap
cd tecap
pip install -e .[dev]
pytest
```

## Quick start

```bash
# Classify reads and emit per-mechanism breakdown
tecap classify \
    --bam sample.bam \
    --gtf gencode.v45.annotation.gtf.gz \
    --polya-sites atlas.clusters.3.0.GRCh38.GENCODE_42.bed.gz \
    --sample S1 \
    --out-dir results/ \
    --threads 8 \
    --platform cdna-pacbio \
    --verbose

# Measure base composition in the 20 nt window downstream of each cleavage site
tecap basecomp \
    --bam sample.bam \
    --gtf gencode.v45.annotation.gtf.gz \
    --polya-sites atlas.clusters.3.0.GRCh38.GENCODE_42.bed.gz \
    --fasta GRCh38.primary_assembly.genome.fa.gz \
    --sample S1 \
    --out-dir results/ \
    --threads 8 \
    --verbose

# Cross-sample comparison
tecap compare \
    --mode classify \
    --inputs results/A_terminal_exon.json,results/B_terminal_exon.json \
    --out-dir results/

tecap compare \
    --mode basecomp \
    --inputs results/A_basecomp.json,results/B_basecomp.json \
    --out-dir results/

# Fetch references (PolyASite 3.0 + matching GENCODE release)
tecap download-atlas \
    --genome GRCh38 \
    --gtf-version 45 \
    --out-dir ref/
```

## Outputs

Per sample (classify):
- `{sample}_terminal_exon.json` — bucket counts, fractions, PAS split, UTR-length stratification, orientation sanity check, read-length medians.
- `{sample}_terminal_exon.png` — 3-panel summary plot.
- `{sample}_mecha_scatter.png` — read length vs TE coverage for MechA-correct reads.
- `{sample}_tecap_mqc.json` — MultiQC custom-content table (auto-detected by the `_mqc.json` suffix).
- `{sample}_per_gene.tsv` (optional, with `--per-gene-table`) — per-gene bucket counts.

Per sample (basecomp):
- `{sample}_basecomp.json` — %A histograms per bucket, medians, >=60% and 30-50% fractions.
- `{sample}_basecomp.png` — 8-panel histogram grid.

Cross-sample:
- `comparison_terminal_exon.png` — grouped bars across samples.
- `comparison_basecomp.png` — per-bucket histogram overlays.

## Citation

If you use `tecap`, please cite the GitHub release DOI (see `CITATION.cff`).

## License

MIT
