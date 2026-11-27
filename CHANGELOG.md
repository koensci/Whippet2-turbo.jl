# Changelog

## v2.0.0-fork.3

Upgraded `whippet-delta.jl` with a multithreaded Beta-Binomial GLM differential splicing engine supporting covariates, Empirical Bayes dispersion shrinkage, and Benjamini–Hochberg FDR (developed using Gemini 4 Neon). Requires zero new Julia dependencies and maintains full `.jls` index and `.psi.gz` compatibility with `v2.0.0-fork.1` and `v2.0.0-fork.2`.

### Index Building Improvements (`whippet-index`)
* **8.2× faster human index building (`233.0 s → 28.5 s`, `12.0×` faster than original `Whippet.jl`, `8.74 GB → 4.88 GB` RSS; `src/index.jl`, `src/graph.jl`, `src/refset.jl`, `src/edges.jl`, `bin/whippet-index.jl`):**
  * Eliminated $O(G^2)$ `xcript *= curgraph.seq` concatenation (~450 GiB of temporary copies on human) by encoding `IntVector{2,UInt8}` directly from `xgraph` (`twobit_enc_graphs`) and parallelizing per-chromosome `SpliceGraph` construction with `Threads.@threads`.
  * Replaced $O(N^2)$ per-gene `seq *= nodeseq` / `reverse_complement(nodeseq) * seq`, `path_to_seq` string round-trips, and per-base slicing in `ExpectedGC` with pre-allocated single-pass `LongSequence` construction and a zero-allocation $O(L)$ sliding-window GC counter.
  * Replaced $O(N \times M)$ `unique_tuple` Cartesian loops and tuple-length JIT specialization in `load_gtf` with per-gene `Vector{CoordInt}` accumulation (`sort!`/`unique!` once per gene) and a zero-allocation GTF tab/attribute scanner.
  * Replaced per-push `sort |> unique` in `add_kmer_edge!` with $O(1)$ monotonic append-dedup. Maintains **100% `.jls` binary compatibility**.

### Batch Cohort Quantification (`whippet-quant`)
* **Single-load `--batch <manifest.csv|tsv>` mode & memory reclamation (`bin/whippet-quant.jl`, `src/bias.jl`):** Quantifies an entire cohort in one Julia session (loading and deserializing the `.jls` index once), supporting both 3-column (`sample,r1,out_prefix`) and 4-column (`sample,r1,r2,out_prefix`) manifests. Added lazy empty sentinels and eager buffer release in `JointBiasCounter` (`primer_adjust!`/`gc_adjust!`), early clearing of `sgq.long`/`multi.map`/`quant.classes`, and inter-sample `GC.gc(true)` + OS malloc trim (`malloc_zone_pressure_relief` on macOS / `malloc_trim` on Linux) to keep peak memory flat across batch samples.

### Bug Fixes
* **Tandem-start (`TS`) path-to-node PSI reversal (`src/io.jl`):** Fixed an upstream bug in `output_utr` where `build_utr_graph` popped `TXST_MOTIF` (`TS`) nodes in descending order (`end_n` down to `start_n`) while `output_utr` iterated `n in st:en` with ascending index `i`, swapping upstream and downstream `TS` PSI/CI/`Inc_Paths` rows. Added a 3-node `TS` unit test in `test/runtests.jl` (all 601,423 assertions pass).

### Differential Splicing Improvements (`whippet-delta`)
* **Beta-Binomial GLM with Cox–Reid Empirical Bayes shrinkage (`src/diff_glm.jl`, `src/Whippet.jl`):**
  * Continuous Beta-Binomial log-likelihood handling fractional read counts (`k = Psi * Total_Reads`) with damped Fisher scoring (IRLS) and light ridge stabilization for complete separation (`0.0` vs `1.0` PSI).
  * 3-pass multithreaded estimation: lockstep `.psi.gz` cohort loading (`read_psi_cohort`), constitutive/boundary event exclusion (`_is_polymorphic_event`), per-event Cox–Reid Adjusted Profile Likelihood (APL) dispersion estimation, event-class stratified Empirical Bayes priors (separate Spliced vs. Tandem UTR `TS`/`TE` priors with Smyth finite-df $\log\chi^2$ median bias correction; disable with `--no-stratify-eb`), and MAP shrunk dispersion (`estimate_eb_theta`).
  * Finite-sample $F(\text{df},\, d_{\text{res}} + d_0)$ likelihood-ratio tests (`edgeR`-style effective prior degrees of freedom) and Benjamini–Hochberg `FDR`, plus covariate-adjusted marginal `Psi_A`, `Psi_B`, `DeltaPsi`, `LogOddsRatio`, `SE`, and `Dispersion`.
* **DESeq2-style covariate formulas & CLI (`bin/whippet-delta.jl`):**
  * Added `--design`, `--file-col` (default `file`), `--formula` (supporting additive terms, continuous numeric covariates, and `*`/`:` interactions), and `--test` (explicit contrasts like `Sample_Type:ewing_vs_msc`, binary factors, continuous covariates, or multi-level omnibus LRTs), with automatic per-event rank/collinearity handling under missingness.
  * Classic `-a`/`-b` mode now uses the new Beta-Binomial GLM by default; added `--no-shrink`, `--fixed-dispersion`, and `--legacy` (to reproduce original Whippet Beta Monte Carlo).
* **Tests (`test/runtests.jl`):** Added 29 unit and end-to-end CLI tests in `@testset "Diff and Delta"` covering Beta-Binomial MLE/LRT, complete separation, BH FDR, formula parsing, covariate adjustment on batch-confounded `.psi.gz` cohorts, continuous covariates, and CLI execution.

---

## v2.0.0-fork.2

Second wave of `whippet-quant` performance optimizations (developed using Gemini 4 Neon). Maintains **100% `.jls` index compatibility** and **byte-for-byte identical uncompressed outputs** with `v2.0.0-fork.1`.

### Performance Improvements (`whippet-quant`)
* **3-stage pipelined paired-end execution (`src/reads.jl`):** Added `PairedChunkBuf` with 3 disjoint ring buffers overlapping Stage 1 (concurrent `R1`/`R2` FASTQ reading), Stage 2 (parallel alignment across workers), and Stage 3 (strictly sequential deterministic counting, bias updates, and optional SAM output).
* **Reusable `SGAlignment` pools & copy-on-insert (`src/align.jl`, `src/paired.jl`, `src/quant.jl`, `src/bias.jl`):**
  * Added `reset_align!`, `_ungapped_align!`, and `ungapped_align_into!` so both single-threaded (`-t 1`) and multi-threaded paired-end alignment reuse pre-allocated `SGAlignment` slots without per-pair heap allocations on uniquely seeded non-branching pairs.
  * Updated `count!` (`sgquant.long`), `Base.push!(::MultiMapping, ...)` (`ambig.map`), and `pushzero!` (`src/bias.jl`) to use single-lookup `Base.ht_keyindex` / `pushto!` with copy-on-insert for new alignment paths.
* **Exact FM-index seed re-check skip (`src/index.jl`, `src/align.jl`):** Added SWAR ambiguity detection (`_has_ambig_dna4`) and a `has_ambig::Bool` flag on the in-memory `FastSpliceGraph` wrapper; when a gene has no ambiguous bases and the seed fits inside `offset_node`, `ungapped_fwd_extend` skips re-checking the exact `seed_length - 1` forward seed bases.
* **Accelerated gzip FASTQ decompression (`src/reads.jl`):** Auto-detects `igzip` (falling back to `pigz`, then `Libz.ZlibInflateInputStream` with a 256 KB buffer) on Unix via a direct POSIX `popen`/`fread` stream (`_PopenReadStream`) + 256 KB `BufferedInputStream`, avoiding Julia libuv `iolock` contention across threads.
* **Parallel `process_events` with concurrent gzip compression (`src/events.jl`):** Partitions genes into ordered blocks, formats and compresses each block concurrently into RFC 1952 concatenated gzip members (`_NoCloseIOBuffer` + `ZlibDeflateOutputStream`), and writes blocks in exact gene order (100% compatible with `gzip`/`gzcat` and `whippet-delta.jl`).

---

## v2.0.0-fork.1

Fork of the `whippet2-julia19-dev` branch of [`timbitz/Whippet.jl`](https://github.com/timbitz/Whippet.jl) (`whippet2` branch). Changes in this release were developed using Gemini 4 Neon.

### Performance Improvements (`whippet-quant`)
* **Concrete in-memory index wrappers (`src/index.jl`, `src/Whippet.jl`, `bin/whippet-quant.jl`):** Added zero-copy `FastGraphLib{K, F}` and `FastSpliceGraph{K}` wrappers with concrete `LongSequence{DNAAlphabet{4}}`, `Edges{K}`, and `FMIndex` fields to eliminate dynamic dispatch in alignment and quantification while maintaining **100% compatibility with existing Julia 1.9 `.jls` index files**.
* **Zero-allocation alignment hot path (`src/record.jl`, `src/fmindex_patch.jl`, `src/sgkmer.jl`, `src/edges.jl`, `src/align.jl`, `src/paired.jl`):**
  * Concretely typed `FASTQRecord.sequence` (`LongSequence{DNAAlphabet{2}}`) with in-place `sequence`/`quality` buffer reuse in `fill!`.
  * Zero-allocation subrange FM-index search (`sa_range_sub` / `sa_range_sub_rc`) and k-mer indexing (`kmer_index_sub`).
  * Immutable `isbits` `SGAlignScore`, scalar register accumulation during extension, shallow alignment copying instead of `deepcopy`, 256-entry `phred_to_prob` lookup table, direct 2-bit vs 4-bit base comparisons, and $O(\log N)$ binary-search edge lookup (`intersect_sorted_gene_range`).
  * Early exit on empty mate seeds and zero-allocation fast path for uniquely seeded read pairs.
* **Multithreaded paired-end & single-end alignment (`src/reads.jl`):** Added `Threads.@spawn` chunked alignment with concurrent double-buffered `R1`/`R2` gzip reading and deterministic count accumulation (use `julia -t auto` or `julia -t N`), and removed manual `GC.gc()` pauses every 100k reads.
* **Bias correction (`src/bias.jl`):** Zero-allocation rolling 2-bit hexamer primer counting and sliding-window GC content counting.
* **Quantification, EM, events & I/O (`src/quant.jl`, `src/events.jl`, `src/io.jl`, `src/graph.jl`):**
  * Skip zero-count nodes/edges in `build_equivalence_classes!` and avoid repeated `collect(iset)` allocations.
  * Track active isoforms and compact converged equivalence classes in `gene_em!`.
  * Canonical fast path in `_process_events`, hoisted `maximum(edges)`, type-stable `_process_spliced`, `Set`-backed deduplication in `reduce_graph`, iterative `rec_tandem_em!`, and interval-indexed `edges_write`.
  * Zero-allocation `write_int` / `write_float` (via `Base.Ryu.writeshortest`) for `.tpm.gz`, `.jnc.gz`, and `.psi.gz` output.
  * Removed obsolete `Base.get(::Tuple, ...)` overwrite in `src/graph.jl` to restore clean Julia 1.9 precompilation.

### Bug Fixes
* **GC bias lookup (`src/bias.jl`):** Fixed `value!(mod::JointBiasMod, bin::Int)` indexing 5'-primer arrays (`mod.fore`/`mod.back`) instead of GC arrays (`mod.gcfore`/`mod.gcback`), and guarded zero denominators.
* **Paired-end GC accumulation (`src/bias.jl`):** Fixed `count!(mod::JointBiasMod, fwd, rev)` overwriting mate 1's GC bins in `mod.temp.gc` with mate 2's, and dynamically grow `mod.temp.gc` for long reads.
* **Short-read & zero-sum safety (`src/bias.jl`, `src/quant.jl`):** Clamped `primer_count!` loops to read length (preventing `BoundsError` on mapped reads < 37 nt with `--biascorrect`), guarded zero sums in `primer_normalize!`/`gc_normalize!`, and fixed `add_to_back!` return value.
* **SRA FASTQ `'N'` encoding (`src/record.jl`):** Sanitize ambiguous bases in-place via a 256-byte lookup table in `Base.fill!`, fixing `Cannot encode 78 to BioSequences.DNAAlphabet{2}()` on SRA FASTQ files with `+SRR...` line-3 headers ([upstream #135](https://github.com/timbitz/Whippet.jl/issues/135)).
* **Tandem-start (`TS`) event indexing (`src/io.jl`, `src/events.jl`):** Fixed `output_utr` comparing `String` `"TS"` to `EdgeMotif` `TXST_MOTIF` (which shifted internal `TS` events with `node > 1` down by 1 node), and updated `next_istandem` to check `== TXEN_MOTIF`.
* **Other fixes:** Guarded zero-total division (`NaN` injection) in `assign_ambig!` (`src/quant.jl`), fixed hardcoded `prop_temp = zeros(500)` limit in `gene_em!` (`src/quant.jl`), fixed supplementary mate-2 `fwd_mate` flag in paired-end `write_sam` (`src/io.jl`), fixed `sgkmer(s)` -> `sgkmer(sub)` in `add_kmer_edge!` (`src/edges.jl`), and fixed duplicate `one` check in `EdgeType` assertion (`src/graph.jl`).

### Compatibility Note
* **Indexes:** Existing Julia 1.9 `.jls` index files work without rebuilding.
* **Differential splicing (`whippet-delta.jl`):** Because the `TS` node-index and `--biascorrect` fixes change output rows/values, all `.psi.gz` files compared together in `whippet-delta.jl` should be quantified with the same version.
