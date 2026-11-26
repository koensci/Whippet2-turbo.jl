# Whippet2-turbo (`Whippet.jl` fork)

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

> **Fork Notice:** **Whippet2-turbo** (`whippet2` branch, `v2.0.0-fork.3`) is a performance-optimized, bug-fixed, and statistically upgraded fork of the `whippet2-julia19-dev` branch of [timbitz/Whippet.jl](https://github.com/timbitz/Whippet.jl). The performance optimizations, bug fixes, and Beta-Binomial GLM differential splicing engine in this fork were developed using Gemini 4 Neon. All original Whippet splice-graph algorithms, architecture, and MIT licensing remain the work of Tim Sterne-Weiler ([Sterne-Weiler et al., *Molecular Cell* 2018](https://doi.org/10.1016/j.molcel.2018.08.018)).

![Whippet Schematic](https://timbitz.github.io/images/Whippet-Fig1.gif)

---

## Highlights of Whippet2-turbo

* **Up to 16× faster `whippet-quant` & 12× faster event PSI calculation** on paired-end FASTQ + `--biascorrect`:
  * Zero-copy concrete in-memory `FastGraphLib` wrappers (**100% binary-compatible with existing Julia 1.9 `.jls` index files** — no index rebuild required).
  * 3-stage pipelined multithreaded paired-end alignment (`julia -t auto`) overlapping concurrent R1/R2 FASTQ decompression, parallel alignment across reusable `SGAlignment` pools, and deterministic sequential counting/bias updates.
  * Automatic `igzip` / `pigz` FASTQ decompression via POSIX `popen` (avoiding libuv lock contention; automatically falls back to `Libz`).
  * Exact FM-index seed re-check skipping, zero-allocation FASTQ/seed/k-mer/bias/output hot paths, and parallel `process_events` with concurrent RFC 1952 gzip compression.
* **New Beta-Binomial GLM `whippet-delta.jl` with covariates & FDR:**
  * Replaces legacy 10,000-draw Monte Carlo Beta comparisons by default with a multithreaded **Beta-Binomial GLM** featuring **Cox–Reid Empirical Bayes dispersion shrinkage** (with finite-df $\log\chi^2$ median bias correction), **finite-sample $F$-test LRT p-values**, and **Benjamini–Hochberg FDR**.
  * Supports **DESeq2-style covariate formulas** (`--design`, `--file-col`, `--formula "~ batch + condition"`, `--test condition:treated_vs_ctrl`) including categorical covariates, continuous numeric covariates, interaction terms (`*`, `:`), and automatic per-event collinearity handling under missingness.
  * Drop-in compatible with classic `-a`/`-b` invocations (and includes `--legacy` to reproduce the original Monte Carlo Beta method).
* **Critical Bug Fixes over upstream `whippet2-julia19-dev`:**
  * Fixed `--biascorrect` GC bias lookup indexing primer arrays (`fore`/`back`) instead of `gcfore`/`gcback`, and fixed paired-end `gc_count!` overwriting mate 1's GC bins.
  * Fixed internal tandem-start (`TS`) node indexing in `output_utr` and `next_istandem` bounds.
  * Fixed SRA FASTQ `'N'` encoding failures (`Cannot encode 78`, [upstream #135](https://github.com/timbitz/Whippet.jl/issues/135)) — raw SRA FASTQ files with `+SRR...` line-3 headers run directly without preprocessing.
  * Guarded short reads (<37 nt) in `--biascorrect`, guarded zero-weight `NaN`s in `assign_ambig!`, fixed supplementary mate-2 SAM flags, and removed hardcoded 500-isoform EM buffer limits.

### Benchmark Summary (50k paired-end 100 bp reads, 58,806-gene human index, `--biascorrect`, warm steady-state)

| Stage | Upstream `whippet2` (1T) | `v2.0.0-fork.3` (1T) | `v2.0.0-fork.3` (8T) | `v2.0.0-fork.3` (16T) | Speedup (16T vs Upstream) |
|---|---:|---:|---:|---:|---:|
| `process_paired_reads!` | 978 ms (412.6 MiB) | 453.5 ms (31.5 MiB) | 92.5 ms (51.9 MiB) | **58.8 ms** (71.3 MiB) | **16.6×** (~850k pairs/s) |
| `process_events` | ~698 ms (67.4 MiB) | 580.1 ms (67.4 MiB) | 86.8 ms (81.9 MiB) | **55.1 ms** (83.5 MiB) | **12.7×** |

See [`CHANGELOG.md`](CHANGELOG.md) for the complete version-by-version history.

---

## Core Features

- **Splice graph representations of transcriptome structure**
  - Build an index for any species from a genome FASTA and GTF annotation.
  - Supplement the index with novel splice sites and exons from aligned RNA-seq BAM files.
  - Dynamic _de novo_ AS event discovery between indexed donor/acceptor splice sites.
- **High-speed spliced read alignment & quantification (`whippet-quant.jl`)**
  - Direct FASTQ-to-splice-graph alignment (single-end or paired-end, read lengths ≤ 255 bp) — no genome BAM alignment required.
  - Multi-mapping gene-family EM assignment, 5'-primer and GC-content bias correction (`--biascorrect`), and optional SAM/circular-RNA output.
  - Node-level PSI, complexity ($K_n$), Shannon entropy, raw junction counts (`.jnc.gz`), and gene/isoform TPMs (`.tpm.gz`).
- **Covariate-aware differential splicing (`whippet-delta.jl`)**
  - Multithreaded Beta-Binomial GLM with Cox–Reid Empirical Bayes dispersion shrinkage, finite-sample $F$-test LRT, BH FDR, and arbitrary covariate design matrices.

---

## Get Started

### 1) Install (Julia 1.9)

Whippet2-turbo runs on **Julia 1.9.x** (available from [julialang.org/downloads/oldreleases](https://julialang.org/downloads/oldreleases/)).

```bash
git clone https://github.com/koensci/Whippet2.jl.git
cd Whippet2.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Verify the installation with the unit test suite (601,418 assertions):
```bash
julia --project=. -t 4 -e 'cd("test"); include("runtests.jl")'
```

*Notes:*
* **Optional speed boost:** Installing `igzip` (`isa-l`) or `pigz` (e.g., `brew install isa-l pigz` on macOS or `apt/conda install isa-l`) is automatically detected by `whippet-quant.jl` for faster gzip FASTQ decompression; otherwise it falls back cleanly to `Libz`.
* Always pass `--project=/path/to/Whippet2.jl` (or `--project=.` from inside the repo) and `-t auto` (or `-t N`) for multithreading.
* For all executables in `bin/`, pass `-h` to view all command-line options and defaults.

---

### 2) Build an Index (`whippet-index.jl`)

> **Compatibility:** Existing Julia 1.9 `.jls` index files built with `whippet2-julia19-dev` work out of the box — no rebuild needed.

#### a) Annotation (GTF) only index
Provide a genome FASTA and an Ensembl/GENCODE-style GTF2.2 file (default hg19/mm10 GTFs are included in `anno/`):

```bash
julia --project=. bin/whippet-index.jl \
  --fasta genome.fa.gz \
  --gtf annotation.gtf.gz \
  -x index/graph
```

*Notes:*
* Whippet uses GTF `exon` lines containing both `gene_id` and `transcript_id` attributes (which must differ from one another), with each transcript's exons in a contiguous block ([GTF2.2](http://mblab.wustl.edu/GTF22.html)). UCSC Table Browser / RefSeq GTFs generally do not meet this specification.
* The default index prefix (`-x`) is `index/graph` (creating `index/graph.jls` and `index/graph.jls.exons.gz`).

#### b) Annotation (GTF) + BAM supplemented index
To incorporate unannotated splice sites and exons discovered in aligned RNA-seq BAM files, first sort and index your BAM with `samtools`:

```bash
samtools sort -o sample.sort.bam sample.bam
samtools index sample.sort.bam

julia --project=. bin/whippet-index.jl \
  --fasta genome.fa.gz \
  --gtf annotation.gtf.gz \
  --bam sample.sort.bam \
  --bam-min-reads 2 \
  -x index/graph_bam
```

*Notes:*
* `--bam` is strand-sensitive; strand-specific alignments are recommended.
* By default, only splice junctions where at least one splice site matches the annotation are added. Pass `--bam-both-novel` to allow fully novel junctions, and tune `--bam-min-reads` (default `1`) to filter low-support cryptic junctions.

---

### 3) Quantify Reads (`whippet-quant.jl`)

#### a) Paired-end reads (recommended with `-t auto` and `--biascorrect`)
```bash
julia --project=. -t auto bin/whippet-quant.jl \
  R1.fastq.gz R2.fastq.gz \
  -x index/graph.jls \
  --biascorrect \
  -o sample_out
```

#### b) Single-end reads
```bash
julia --project=. -t auto bin/whippet-quant.jl \
  sample.fastq.gz \
  -x index/graph.jls \
  --biascorrect \
  -o sample_out
```

#### c) Additional options
* **SAM output:** Add `--sam` to stream alignments to `stdout` (e.g., `... --sam | samtools view -bS - > sample.bam`).
* **Circular RNA back-splicing:** Add `--circ` to quantify back-spliced (`BS`) junctions.
* **Pooled / piped gzip streams without `.gz` extension:** Add `--force-gz`:
  ```bash
  julia --project=. -t auto bin/whippet-quant.jl <( cat lane_{1,2}.fastq.gz ) --force-gz -x index/graph.jls -o pooled_out
  ```
* **Raw SRA FASTQ files:** Files downloaded via `fasterq-dump`/`fastq-dump` with `+SRR...` line-3 headers and `'N'` bases are handled natively without preprocessing.

---

### 4) Differential Splicing (`whippet-delta.jl`)

`whippet-delta.jl` fits a **Beta-Binomial GLM** across samples using each event's inclusion count ($k_{ij} = \Psi_{ij} N_{ij}$) and total read depth ($N_{ij}$), with **Cox–Reid Empirical Bayes dispersion shrinkage**, **finite-sample $F$-test LRT p-values**, and **Benjamini–Hochberg FDR**. See [`docs/GLM_METHODS.md`](docs/GLM_METHODS.md) for the full statistical formulation, equations, and example manuscript methods text.

#### a) Simple 2-group comparison (`-a` vs `-b`)
Pass glob patterns or comma-separated file lists for Set A and Set B (`DeltaPsi = Psi_A - Psi_B`):

```bash
julia --project=. -t auto bin/whippet-delta.jl \
  -a treat_ -b ctrl_ \
  -r 10 -s 2 \
  -o treat_vs_ctrl

# Or explicit comma-separated file lists:
julia --project=. -t auto bin/whippet-delta.jl \
  -a trt1.psi.gz,trt2.psi.gz,trt3.psi.gz \
  -b ctl1.psi.gz,ctl2.psi.gz,ctl3.psi.gz \
  -o treat_vs_ctrl
```

#### b) Covariate-adjusted GLM (`--design`, `--formula`, `--test`)
Provide a CSV or TSV metadata table containing a file-path column (default `file`, or specify any column name via `--file-col`) plus any categorical or numeric covariates:

Example `metadata.csv`:
```csv
sample,whippet_output,Sample_Type,Stranded,batch,age
s1,/path/to/s1.psi.gz,ewing,stranded,B1,16
s2,/path/to/s2.psi.gz,ewing,unstranded,B2,21
s3,/path/to/s3.psi.gz,msc,stranded,B1,28
s4,/path/to/s4.psi.gz,msc,unstranded,B2,34
```

Run a covariate-adjusted contrast (`Sample_Type:ewing_vs_msc` sets `msc` as reference and computes `DeltaPsi = Psi(ewing) - Psi(msc)`):

```bash
julia --project=. -t auto bin/whippet-delta.jl \
  --design metadata.csv \
  --file-col whippet_output \
  --formula "~ Stranded + Sample_Type" \
  --test Sample_Type:ewing_vs_msc \
  -r 10 -s 5 \
  -o ewing_vs_msc_glm
```

**Supported formula & test syntax:**
* **Additive covariates:** `--formula "~ batch + donor + condition"` (categorical factors and numeric continuous covariates are auto-detected).
* **Interactions:** `--formula "~ batch * condition"` (expands to `batch + condition + batch:condition`) or explicit `--formula "~ batch + condition + batch:condition"`.
* **Explicit 2-level contrast:** `--test condition:treated_vs_ctrl` (works on both binary and 3+ level factors, reporting covariate-adjusted marginal `Psi_A` and `Psi_B`).
* **Continuous predictor:** `--formula "~ batch + age" --test age` (reports average marginal slope `DeltaPsi = dPsi/dx`, `LogOddsRatio`, `SE`, `Pvalue`, `FDR`).
* **Multi-level omnibus LRT:** `--test condition` on a factor with $K > 2$ levels runs a $(K-1)$-df LRT.

#### c) Key `whippet-delta.jl` options
| Flag | Default | Description |
|---|---|---|
| `-r, --min-reads` | `5` | Minimum `Total_Reads` in a sample for that sample to be included for an event. |
| `-s, --min-samples` | `1` | Minimum valid samples required (per group in 2-group/contrast mode; total in general GLM mode). |
| `--file-col` | `file` | Column name in `--design` containing `.psi.gz` paths. |
| `--no-shrink` | off | Disable Empirical Bayes dispersion shrinkage (use per-event Cox–Reid dispersion MLE). |
| `--fixed-dispersion` | off | Fix Beta-Binomial overdispersion $\phi \in (0, 1)$ (recommended for unreplicated 1 vs 1 exploratory runs, e.g., `--fixed-dispersion 0.05`). |
| `--legacy` | off | Use original Whippet 10,000-draw Beta Monte Carlo method (`-a`/`-b` mode only). |

*Note on collinearity:* Ensure tested variables are not completely confounded with covariates (e.g., if every study or sequencer belongs exclusively to one condition, do not include it in `--formula`). Per-event missingness that makes an untested covariate level disappear is handled automatically via pivoted QR rank reduction.

---

## Output Formats

### 1) `.gene.tpm.gz` and `.isoform.tpm.gz` (`whippet-quant.jl`)

| Gene/Isoform | TpM | Read_Counts |
|---|---:|---:|
| NFIA | 2897.11 | 24657.0 |

### 2) `.psi.gz` (`whippet-quant.jl`)

Quantified splice-graph nodes fall into these event types:

| Type | Interpretation |
|---|---|
| `CE` | Core exon (may be bounded by one or more alternative `AA`/`AD` nodes) |
| `AA` | Alternative Acceptor splice site |
| `AD` | Alternative Donor splice site |
| `RI` | Retained intron |
| `TS` | Tandem transcription start site |
| `TE` | Tandem alternative polyadenylation site |
| `AF` | Alternative First exon |
| `AL` | Alternative Last exon |
| `BS` | Circular back-splicing (output only when `--circ` is given) |

![Whippet AS Events](https://timbitz.github.io/images/Whippet-Event-Types.gif)

| Gene | Node | Coord | Strand | Type | Psi | CI_Width | CI_Lo,Hi | Total_Reads | Complexity | Entropy | Inc_Paths | Exc_Paths | Edges |
|---|---:|---|---|---|---|---|---|---|---|---|---|---|---|
| NFIA | 2 | chr1:61547534-61547719 | + | AF | 0.782 | 0.191 | 0.669,0.86 | 49.0 | K1 | 0.756 | 2-4-5:0.782 | 1-5:0.218 | ... |
| NFIA | 4 | chr1:61548433-61548490 | + | CE | 0.8329 | 0.069 | 0.795,0.864 | 318.0 | K2 | 1.25 | 2-4-5:0.3342,3-4-5:0.4987 | 1-5:0.1671 | ... |

Whippet quantifies splicing at the **node** level (rather than enumerating all full-exon combinations), which scales cleanly to arbitrarily complex events ($K_n = \lceil\log_2(\text{paths})\rceil$) and computes Shannon entropy over path probabilities. Use the companion `*.exons.gz` file produced by `whippet-index.jl` when mapping nodes back to full exon coordinates.

### 3) `.jnc.gz` (`whippet-quant.jl`)

| CHROM | DONOR_COORD | ACCEPTOR_COORD | GENEID:NODES:TYPE | READCOUNT | STRAND |
|---|---:|---:|---|---:|---|
| chr1 | 150600068 | 150601890 | ENSG00000143420.17_1:2-3:CON_ANNO | 113.04 | - |
| chr1 | 161185160 | 161188031 | ENSG00000158869.10_1:1-3:ALT_UNIQ | 3.11 | + |

(`CON_ANNO` = annotated constitutive, `ALT_ANNO` = annotated alternative, `ALT_UNIQ` = unannotated alternative.)

### 4) `.diff.gz` (`whippet-delta.jl`)

#### Default Beta-Binomial GLM — 2-group or explicit 2-level contrast mode:
| Column | Meaning |
|---|---|
| `Gene`, `Node`, `Coord`, `Strand`, `Type` | Splice-graph event identifier |
| `Complexity`, `Entropy` | Maximum complexity ($K_n$) and Shannon entropy across valid samples |
| `SampleCnt_A`, `SampleCnt_B` | Number of valid samples (`Total_Reads >= min-reads`) in Group A and Group B |
| `MeanReads_A`, `MeanReads_B` | Mean read depth across valid samples in Group A and Group B |
| `Psi_A`, `Psi_B` | Covariate-adjusted marginal mean PSI for Group A and Group B |
| `DeltaPsi` | `Psi_A - Psi_B` |
| `LogOddsRatio`, `SE` | Fitted logistic regression coefficient ($\hat\beta$) and Fisher standard error |
| `Dispersion` | Empirical Bayes shrunk Beta-Binomial overdispersion $\tilde\phi \in (0, 1)$ |
| `LRT`, `Pvalue`, `FDR` | Likelihood-ratio statistic, finite-sample $F$-test p-value, and Benjamini–Hochberg adjusted FDR |

#### General GLM mode (continuous covariate or multi-level factor omnibus test):
Outputs `Gene`, `Node`, `Coord`, `Strand`, `Type`, `Complexity`, `Entropy`, `SampleCnt`, `MeanReads`, `MeanPsi`, `DeltaPsi` (average marginal slope $d\Psi/dx$ for continuous predictors, `NA` for multi-df omnibus tests), `LogOddsRatio`, `SE`, `Dispersion`, `LRT`, `DF`, `Pvalue`, `FDR`.

#### Legacy mode (`--legacy`):
Outputs the original Whippet format (`... SampleCnt_A SampleCnt_B MeanReads_A MeanReads_B Psi_A Psi_B DeltaPsi Probability`).

---

## Important Compatibility Note

* **`whippet-quant` outputs** are byte-for-byte identical across `v2.0.0-fork.1`–`v2.0.0-fork.3` and use the exact same `.jls` index format as upstream `whippet2-julia19-dev`.
* **Do not mix `.psi.gz` files** generated by Whippet2-turbo with `.psi.gz` files generated by unpatched upstream Whippet in `whippet-delta.jl`, because the `--biascorrect` and internal `TS` node-index bug fixes intentionally change affected event rows/values.

---

## Index Building Strategies & Visualization

* **High-quality annotations:** Prefer curated annotations (e.g., GENCODE TSL-1/2 for human) and strand-specific BAMs when supplementing indexes to avoid cross-gene artifacts in overlapping loci.
* **Short reads (~36 nt):** Consider lowering the index k-mer size (e.g., `6` instead of default `9`).
* **Sashimi plots:** Use `whippet-quant.jl --sam` and convert to BAM with `samtools` to inspect alignments and sashimi plots in IGV.

---

## Citation

If you use Whippet2-turbo in published work, please cite the original Whippet paper and reference this fork/release tag in your Methods:

> Sterne-Weiler T, Weatheritt RJ, Best AJ, Ha KCH, Blencowe BJ. **Efficient and Accurate Quantitative Profiling of Alternative Splicing Patterns of Any Complexity on a Laptop.** *Molecular Cell*. 2018;72(1):187–200.e6. https://doi.org/10.1016/j.molcel.2018.08.018
