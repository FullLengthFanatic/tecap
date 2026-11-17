# tecap: a whitepaper

**3' terminal exon capture diagnostics for long-read single-cell RNA-seq.**

Version 0.4.0. MIT licensed. [github.com/FullLengthFanatic/tecap](https://github.com/FullLengthFanatic/tecap) · [doi:10.5281/zenodo.19762736](https://doi.org/10.5281/zenodo.19762736)

---

## Abstract

`tecap` reads a long-read RNA-seq alignment file and asks a single question of every transcript: where did its 3' end actually come from? Each multi-exon read is sorted into one of nine mechanism buckets based on where its 3' end lands relative to the gene's terminal exon, its 3' UTR, and a reference atlas of real polyadenylation sites. It then counts the adenines in the 20 nt of genome just downstream of each read's site and labels the site templated (8 or more), non-templated (4 or fewer) or ambiguous (5 to 7): the bucket says where the site is, the label says what kind of site it is. It also reports two descriptive bands of the adenine share, 30% to 50% A and 60% A or more. On real data single-cell library preps are enriched in the 30% to 50% band and bulk Iso-Seq in the 60% or more band, a clean split whose biochemical cause is not yet known. tecap is a diagnostic, not a filter: it tells you *why* a protocol's reads end where they do, not whether they pass.

---

## 1. The problem

Almost every RNA-seq library prep that targets mRNA uses **oligo-dT**: a short stretch of DNA T's that anneals to the poly(A) tail at the 3' end of a transcript and primes reverse transcription. When it works, the priming site is the real cleavage-and-polyadenylation site, and the read's 3' end marks where the cell chose to end that transcript.

It does not always work. Genomes are full of A-rich stretches that have nothing to do with polyadenylation: A-runs inside 3' UTRs, inside coding sequence, inside introns. Oligo-dT will anneal to those too. This is **internal priming**, and it produces reads whose 3' end is an artifact: a fake polyA site, somewhere in the middle of a gene. In short-read 3' counting this mostly inflates counts; in long-read single-cell RNA-seq it does something more visible. The read looks like a truncated or mislocated transcript, and the protocol appears to "capture" the true 3' end less often than it should.

That last sentence is the concrete observation tecap was built to dissect. Matched human retinal organoid samples were processed in parallel, one on 10x Chromium v4 and one on an in-house FLASH-seq variant adapted to the BD Rhapsody platform, both sequenced with PacBio Kinnex. Both chemistries use oligo-dT capture and 16 PCR cycles. Despite that, 10x captured the full 3' terminal exon far less often than the BD prep. The gap was unexpected, and "10x is just worse" is not a mechanism. tecap exists to replace that hand-wave with a per-read accounting of where the missing reads went.

### A one-paragraph biology primer

The **terminal exon (TE)** is the last exon of a transcript. It carries the **3' UTR** (the untranslated region after the stop codon) and, at its far end, the real polyA site. A read that captures the true mRNA 3' end will have its 3' end sitting in the TE, at or near an annotated polyA cluster, having covered most of the TE on the way there. A read born of internal priming will have its 3' end somewhere else: short of the polyA site, in the coding sequence, on an upstream exon, in an intron. tecap's whole logic is geometric. Tell it where the 3' end landed and it tells you which of those stories the read is telling.

---

## 2. What tecap does

The core idea is one sentence: **in an oligo-dT prep, a read's priming site is where its terminal adenine run starts.** tecap calls that point the site: aligned and clipped adenines are counted together (10 or more, one other base allowed), because the aligner can align adenines encoded in the genome and carry the alignment end past them. A read without such a run, and every read under `--site-anchor alignment_end`, uses the alignment end instead. Below, a read's 3' end means its site. Its position relative to the gene model is direct evidence about where the read was primed.

tecap turns that into two measurements.

**Classification (the `classify` step).** Every multi-exon read is assigned to exactly one of nine buckets by where its 3' end falls relative to the TE, the TE's CDS/UTR split, and the nearest annotated PolyASite cluster. The buckets range from clean success (the read covers most of the TE) through real-but-truncated, real alternative polyadenylation, and the various flavors of internal priming.

**Base composition (the `basecomp` step).** For every spliced read in the eight basecomp buckets (MechA-internalUTR, MechA-noCDS and MechC are left out), tecap fetches the 20 nucleotides of *reference genome* immediately 3' of the site (in mRNA orientation) and counts the adenines. With 8 or more the site is templated (internal priming), with 4 or fewer non-templated (priming on the poly(A) tail), with 5 to 7 ambiguous; a read whose site fell back to the alignment end gets no label. Group those by bucket and you get a fingerprint of the sequence oligo-dT was annealing to. A real polyA cleavage site is not A-rich downstream; an internal-priming site usually is.

The framing matters and is deliberate. tecap does not call reads "good" or "failed"; its labels say what kind of site a read has, not whether the read passes. A read with a templated site is not a failed read; it is a read that was redirected to a fake polyA site, and which fake site tells you about the chemistry. The output is a mechanism breakdown, read as a diagnostic.

### The gene model and where reads land

For a plus-strand gene, the picture is a left-to-right track from 5' to 3'. The read's 3' end (its priming site) can land in any zone:

```
 5'                                                                      3'   (+ strand)
   [ exon1 ]==intron==[ exon2 ]==intron==[========= terminal exon =========]----->
       :                  :        :        |                  |          |  :
       :                  :        :        |<----- CDS ------>|<- 3'UTR ->|  :
       :                  :        :        :                  :     ^polyA   :
       :                  :        :        :                  :      site    :
   ----+--------+---------+--------+--------+------------------+------+-------+------>
   3'-end zone:                                                              
     intron / gene flank ............................................. MechB-aspecific
     on an upstream exon (no polyA) .................................. MechB-exon
     upstream exon, AT a polyA cluster ............................... MechB-APA  (biological)
     inside TE, in the CDS portion ................................... IP-TE-CDS
     inside TE 3'UTR, AT a polyA cluster ............................. MechA-correct
                                            (and Captured if it covered >=50% of the TE)
     inside TE 3'UTR, NOT at a polyA cluster ......................... MechA-internalUTR
     downstream, past the TE end ..................................... MechC
```

Minus-strand genes are the mirror image. tecap takes the strand from the gene, not from the read, so it looks for the 3' end at the read's left end on a minus-strand gene and at its right end on a plus-strand gene.

### The decision tree

The classifier is deterministic and runs in a fixed order (`src/tecap/classify.py`):

1. If the read covers **>= 50%** of the terminal exon, it is **Captured**. Done.
2. Otherwise, if the 3' end is **upstream** of the TE:
   - at an annotated polyA cluster, it is **MechB-APA** (a real alternative-polyadenylation isoform);
   - else on an upstream exon, it is **MechB-exon**;
   - else (intron or gene flank), it is **MechB-aspecific**.
3. Otherwise, if the 3' end is **downstream** of the TE end, it is **MechC**.
4. Otherwise the 3' end is inside the TE:
   - if the gene is non-coding, it is **MechA-noCDS**;
   - else if it is in the 3' UTR (past the CDS end): at a polyA cluster it is **MechA-correct**, otherwise **MechA-internalUTR**;
   - else it is in the CDS, so **IP-TE-CDS**.

### The nine buckets

The text below is the tool's own definition for each bucket (`tecap explain`).

| Bucket | What it means | Why it matters |
|---|---|---|
| **Captured** | Read 3' end lands in the terminal exon and the read covers >=50% of it. | Successful full-length capture of the mRNA 3' end; the goal of any 3'-end protocol. |
| **MechA-correct** | Read 3' end is in the TE 3' UTR within +-25 bp of an annotated PolyASite cluster, but the read covers <50% of the TE. | Truncated transcript that nonetheless terminates at a real polyA site; common with degraded input or short-fragment library prep. |
| **MechA-internalUTR** | Read 3' end is in the TE 3' UTR but not within +-25 bp of any annotated polyA cluster. | Candidate internal priming site in the UTR; no annotated polyA cluster supports it. |
| **IP-TE-CDS** | Read 3' end is inside the terminal exon's CDS portion (i.e. before the UTR begins). | Candidate internal priming site; real polyA sites are rare in coding sequence. |
| **MechA-noCDS** | Read 3' end is inside the TE of a non-coding gene (no CDS to delimit a UTR). | Coding/UTR split does not apply; reported separately so the coding-gene buckets stay clean. |
| **MechB-APA** | Read 3' end is upstream of the terminal exon but at an annotated polyA cluster on an upstream exon. | Alternative polyadenylation isoform, supported by an annotated polyA cluster. |
| **MechB-exon** | Read 3' end is on an upstream exon, with no annotated polyA cluster nearby. | Candidate internal priming site on an upstream exon; no annotated polyA cluster supports it. |
| **MechB-aspecific** | Read 3' end is upstream of the TE in an intron or gene flank. | Pre-mRNA priming or off-target alignment. |
| **MechC** | Read 3' end is downstream of the terminal exon end. | Read-through, unannotated 3' UTR extension, or alignment artifact past the annotated TE. |

For the polyA-site buckets (MechA-correct, MechB-APA), tecap further splits reads by whether the supporting atlas cluster carries a canonical AAUAAA-like hexamer: **PAS+** (it does, high confidence) versus **PAS-** (it does not). A read at a PAS- cluster is weaker evidence for real cleavage and a candidate mispriming artifact.

---

## 3. How it works

### Inputs

- **BAM** of long-read alignments, coordinate-sorted and indexed. tecap looks only at multi-exon reads (it requires an `N` operation in the CIGAR, i.e. a spliced alignment) and skips unmapped, secondary, and supplementary records. Designed for PacBio Iso-Seq / Kinnex and Oxford Nanopore cDNA.
- **GTF** gene annotation (GENCODE or Ensembl, gzip fine). tecap parses exon and CDS features and derives, per gene, the terminal exon coordinates, the CDS end, the 3' UTR length, and an interval tree of the upstream exons.
- **PolyASite atlas** as a BED file. By default tecap keeps cluster types `TE,DS,DI,EX` and drops the others (`IN,UI,AI,AE,AL,IG`), which carry internal-priming artifacts inherited from the datasets used to build the atlas. Each cluster is flagged PAS+/PAS- from its hexamer field.
- **Indexed reference FASTA** (only for `basecomp`), to read the downstream window.

References are fetched and cached automatically on first run when you pass `--genome` and `--gtf-version`, with a SHA256 check; PolyASite v3.0 for human GRCh38, v2.0 for mouse GRCm38.

Direct-RNA sequencing is refused at the command line on purpose. There is no reverse transcription and no oligo-dT step in dRNA, so there is no priming artifact to diagnose.

### The measurements and their knobs

- **Terminal exon coverage**: overlap of the read with the TE divided by TE length. The default **0.5** threshold is what separates Captured from the truncated buckets.
- **PolyA matching**: a 3' end counts as "at a polyA site" if it falls within **+-25 bp** of an atlas cluster on the same strand.
- **Downstream base composition**: the adenines in the **20 nt** reference window 3' of the site, counted as an integer. Per bucket tecap reports how many reads with an observed site are templated (8 or more), non-templated (4 or fewer) or ambiguous (5 to 7), and the adenine share as a histogram and three bands. `frac_ge_60`, `frac_30_50` and `frac_ge_40` are descriptive bands of the adenine share in the window (60% or more, 30% to under 50%, 40% or more). They do not name a priming mechanism. The label is in `frac_templated`, `frac_nontemplated` and `frac_ambiguous`.
- **UTR-length stratification**: capture and MechA-correct rates are broken out by 3' UTR length bins (0-500, 500-1k, 1k-2k, 2k-5k, >5k bp), because mispriming bias concentrates in long-UTR genes.
- **Orientation sanity check**: the fraction of reads whose alignment strand matches the gene strand (`orient_match_frac`), a check that the input is properly stranded.

### Outputs

Per sample, `classify` writes a JSON (bucket counts and fractions, PAS split, UTR-bin stratification, read-length medians, orientation fraction; schema version 0.1), a three-panel summary PNG, a MechA-correct scatter PNG, and a MultiQC custom-content JSON. An optional per-gene TSV gives bucket counts per gene. `basecomp` writes a per-bucket %A histogram JSON and an eight-panel histogram PNG.

`tecap report` renders a **single self-contained HTML file** with embedded plots, an executive summary, the mechanism legend, per-bucket tables, the PAS split, UTR-length stratification, and a glossary. No JavaScript, no external CSS, no CDN. It runs for one sample or side-by-side across many. `tecap compare` produces just the cross-sample comparison plots, and `tecap explain` prints the bucket definitions at the terminal.

### Engineering

The whole tool is pure Python (`pysam`, `intervaltree`, `numpy`, `matplotlib`); no samtools or bedtools at runtime. Parallelism is by contig: the parent process builds the gene and polyA indices, places them in module globals, then forks a worker pool, so workers inherit the large reference structures by copy-on-write instead of pickling them. Contigs with no mapped reads are skipped before any worker is dispatched. Classification is fully deterministic, with no random sampling anywhere.

---

## 4. What the output reveals: a worked example

### Single-cell vs single-cell (10x vs BD FLASH-seq)

On the matched retinal organoid samples (52.9M multi-exon reads for 10x, 62.3M for BD), the headline gap is real and survives atlas cleanup:

| Category | 10x | BD |
|---|---|---|
| Captured (>=50% TE) | 31.4% | 46.7% |
| MechA-correct | 16.8% | 19.4% |
| MechB-APA | 18.8% | 9.4% |

Three things came out of the per-read accounting that "10x is worse" would never have told you:

- **The PAS+/PAS- split exposes a systematic bias.** Among MechB-APA reads, only 33.4% of 10x reads sit at a hexamer-carrying cluster versus 47.3% for BD. 10x lands disproportionately on PAS- clusters, the ones most likely to be atlas-level internal-priming artifacts.
- **Read length explains about half the gap.** Median aligned length is shorter for 10x (Captured 1036 bp vs 1404 bp). Stratified by UTR length, the capture gap widens monotonically toward long UTRs, exactly what a length-limited read would do. But the gap at short UTRs is too small for length to be the whole story.
- **Windows of 60% A or more are rare.** The fraction of reads with 60% A or more downstream is below 8.5% in every bucket of both samples. What separates the two is the 30% to 50% A band: 10x is shifted warmer there across nearly every bucket, most strongly in MechB-aspecific (+10.5 percentage points, 35.8% vs 25.3%).

### Single-cell vs bulk, across six datasets

Extending the panel to six human datasets on the same backend turns that shift into a clean split. The diagnostic bucket is MechB-aspecific:

| Sample | Prep | n | frac[30,50] | frac>=60 |
|---|---|---:|---:|---:|
| 10x_FL_v02_full | single-cell, GEM droplet | 6,149,697 | 0.358 | 0.029 |
| BD46_FS_SEQ | single-cell, microwell (BD) | 4,630,753 | 0.253 | 0.016 |
| MDA_argentag_kinnex12x | single-cell, nanowell | 2,176,713 | 0.248 | 0.066 |
| PBMC_FS_ONT | single-cell, plate (ONT) | 2,370,853 | 0.298 | 0.072 |
| kinnex_cerebellum | bulk Iso-Seq | 127,951 | 0.164 | 0.468 |
| kinnex_heart | bulk Iso-Seq | 131,125 | 0.160 | 0.475 |

Two families fall out, robust across capture architectures:

1. **Single-cell preps** (10x droplet, BD microwell, ArgenTag nanowell, plate FLASH-seq): the 30% to 50% A band dominates (`frac[30,50]` 0.25 to 0.36), windows of 60% A or more are rare (`frac>=60` at most 0.07).
2. **Bulk Iso-Seq**: the 60% A or more band dominates (`frac>=60` around 0.47), the 30% to 50% A band sits near the genomic baseline.

The separation is roughly 7x on `frac>=60` and 2x on `frac[30,50]`, and it holds across droplet, microwell, nanowell, and plate capture, and across both bead-immobilized and in-solution reverse transcription.

### The honest conclusion

What separates the two families biochemically is not established. Several plausible mechanisms were tested against the panel and ruled out:

- **Saturating local oligo-dT concentration.** Falsified: bulk Iso-Seq runs the same 1.2 uM final oligo-dT as FS-ONT and BD.
- **Template-switching reverse transcription.** Falsified: PacBio Targeted Iso-Seq also template-switches.
- **PCR primer asymmetry.** Falsified: FS-ONT and 10x both use two different amplification primers.
- **RT polymerase identity.** Implausible at this magnitude.

Lysis chemistry, PCR cycle count, and read orientation were matched and are not the axis. Until the underlying biochemistry is identified, the split between the 30% to 50% A band and the 60% A or more band is reported as an empirical signature that distinguishes single-cell preps from bulk Iso-Seq, with no claim about cause. The tool's value is that it surfaces and quantifies this axis cleanly; the biology is a separate, open question.

---

## 5. Strengths

- **It answers a mechanism question, not a pass/fail one.** A 15-point capture gap becomes a labeled accounting of where the reads went and how they were primed. That is actionable in a way a single QC number is not.
- **Geometry is precise and strand-aware.** Coordinate conventions (1-based GTF to 0-based half-open intervals), terminal-exon selection, and the upstream/downstream tests are handled consistently for both strands and are unit-tested across all nine buckets.
- **It fights annotation noise instead of trusting it.** Filtering the atlas to real polyA cluster types and splitting every polyA-site call into PAS+/PAS- both push back on the internal-priming artifacts baked into reference atlases. The 10x MechB-APA inflation (much of it PAS-) is a direct example.
- **Base composition is read off the reference, not the read.** That makes the A-content measurement independent of basecalling error near the cleavage site, which matters for ONT data.
- **The reporting is built for people who did not run the tool.** Self-contained HTML, an intro paragraph, a reading guide, and a glossary mean a reviewer or collaborator can open one file and understand it. Definitions live in one place (`constants.py`) and feed the README, the reports, the plot captions, and `tecap explain` identically.
- **It is light to deploy and run.** Pure Python, no external CLI tools, contig-level parallelism with copy-on-write references, and a low memory footprint.
- **It is scientifically honest in its own record.** An earlier version attributed the 30% to 50% A split to a specific oligo-dT mechanism; when bench data falsified that, the claim was retracted in the changelog and the published release notes were left intact rather than quietly rewritten.

---

## 6. Limitations

These are limitations of the method and its scope, separate from the engineering roadmap.

- **The headline finding has no mechanism.** The single-cell vs bulk split is empirical. tecap measures it; it does not explain it, and it should not be cited as if it did.
- **Base composition is reference-based.** It reads genomic A-content, not the actual nucleotides the oligo-dT annealed to in a given molecule. Polymorphic adenine runs or errors near the cleavage point can shift the real target. Over millions of reads this averages out; for individual loci it does not.
- **The 20 nt window is short.** A 30 or 40 nt window might resolve the 30% to 50% A signature more cleanly. The window size is configurable; the 30/50/60% A band edges are not, currently, and the labels are defined for the 20 nt window only (with any other `--window` the label fields are null).
- **The 50% coverage threshold is arbitrary.** It sets the Captured vs MechA-correct boundary. Bucket *ordering* across samples is stable to the choice, but absolute percentages move with it.
- **The polyA atlas is imperfect.** PolyASite inherits internal-priming artifacts from its source data. The PAS+/PAS- split mitigates this but does not eliminate it: a genuine non-canonical polyA event can be PAS-.
- **Stranded input is assumed, and not enforced.** tecap uses the gene strand to define "downstream" and to look for the read's terminal adenine run. An unstranded BAM dilutes the signal toward the genomic mean. The `PBMC_FS_ONT` sample (orientation match 0.495, meaning random strand because pychopper was not run upstream) is a worked example: its basecomp numbers are conservative bounds, not point estimates. tecap reports the orientation match fraction, and in the development version after 0.4.0 both `classify` and `basecomp` warn when it drops below 0.80, the signature of unstranded input. The guard surfaces the problem; it does not fix it. tecap still writes the diluted metrics rather than correcting or refusing the run, so orienting reads with pychopper / cdna_classifier upstream remains the real fix.
- **It depends on annotation.** Genes or 3' ends absent from the GTF, or transcripts with unannotated 3' UTR extensions, are mis-bucketed (often into MechC). Novel-isoform discovery is out of scope.
- **Only multi-exon reads are analyzed.** Single-exon transcripts and intronless genes are not classified.
- **Direct-RNA is unsupported by design.** No RT, no oligo-dT, nothing to diagnose.
- **One coverage gap is open.** As of the v0.4 work there was no clean public human plate-based Smart-seq + PacBio HiFi dataset to test whether the 30% to 50% A signature tracks reaction-volume scale or chemistry lineage. That experiment is still missing.

---

## 7. Availability and citation

```bash
pip install git+https://github.com/FullLengthFanatic/tecap@v0.4.0
```

A typical run, with references auto-fetched and cached on first use:

```bash
tecap classify --bam sample.bam --genome GRCh38 --gtf-version 45 \
    --sample S1 --out-dir results/ --threads 8 --platform cdna-pacbio

tecap basecomp --bam sample.bam --genome GRCh38 --gtf-version 45 \
    --fasta GRCh38.primary_assembly.genome.fa.gz \
    --sample S1 --out-dir results/ --threads 8

tecap report --classify-json results/S1_terminal_exon.json \
    --basecomp-json results/S1_basecomp.json --out-html results/S1_report.html
```

tecap is MIT licensed. If you use it, cite the GitHub release DOI ([10.5281/zenodo.19762736](https://doi.org/10.5281/zenodo.19762736), see `CITATION.cff`).

---

## Glossary

| Term | Definition |
|---|---|
| **TE** | Terminal exon. The last exon of a transcript, where the 3' UTR and the canonical polyA site live. tecap classifies reads by where their 3' end lands relative to the TE. |
| **UTR** | Untranslated region. The 3' UTR is the part of the terminal exon downstream of the stop codon, between the CDS end and the polyA site. Most real polyA cleavage happens here. |
| **CDS** | Coding sequence. The protein-coding portion of a transcript. A 3' end inside the TE's CDS portion is flagged IP-TE-CDS, a strong internal-priming signal. |
| **polyA** | Polyadenylation. The 3' processing step that cleaves the pre-mRNA and adds a poly(A) tail. tecap uses the PolyASite atlas to know where real cleavage clusters are. |
| **PAS** | Polyadenylation signal. The canonical AAUAAA-like hexamer ~10-30 nt upstream of a real cleavage site. tecap splits polyA-site reads into PAS+ (hexamer present) and PAS- (absent). |
| **APA** | Alternative polyadenylation. mRNA isoforms that differ in their 3' end through use of upstream polyA sites. MechB-APA captures these biological events and separates them from artifacts. |
| **oligo-dT** | A short DNA oligo of consecutive T's used to capture poly(A)-tailed mRNA. Internal priming happens when it anneals to A-rich genomic stretches instead of the real poly(A) tail. |
