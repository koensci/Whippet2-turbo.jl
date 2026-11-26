#!/usr/bin/env julia
# Tim Sterne-Weiler 2015

using Pkg

const dir = abspath( splitdir(@__FILE__)[1] )
const ver = readline(open(dir * "/VERSION"))

const minversion = v"1.7"

start = time_ns()
println( stderr, "Whippet $ver loading... " )

Pkg.activate(dir * "/..")

using Whippet
using Serialization
using ArgParse

function parse_cmd()
   s = ArgParseSettings()

   @add_arg_table! s begin
      "filename.fastq[.gz]--OR--pre-aligned.bam"
         arg_type = String
         required = false
      "paired_mate.fastq[.gz]"
         arg_type = String
         required = false
      "--batch"
         help     = "Batch manifest CSV/TSV (sample,r1,r2,out_prefix or sample,r1,out_prefix); loads index once across all samples"
         arg_type = String
         default  = ""
      "--index", "-x"
         help     = "Output prefix for saving index 'dir/prefix' (default Whippet/index/graph[.jls])"
         arg_type = String
         default  = fixpath( "$dir/../index/graph" )
      "--out", "-o"
         help     = "Where should the gzipped output go 'dir/prefix'?"
         arg_type = String
         default  = fixpath( "./output" )
      "--biascorrect"
         help     = "Apply fragment GC-content and 5' Sequence bias correction methods for more stable PSI values"
         action   = :store_true

      "--aberrant-min-overlap"
         help     = "Minimum number of overlap in the first or last node of a path (.bam input only)"
         arg_type = Int
         default  = 9
      "--aberrant-min-reads"
         help     = "Minimum number of reads supporting an aberrant splicing event (.bam input only)"
         arg_type = Int
         default  = 2
      "--aberrant-max-xce"
         help     = "Maximum length of a cryptic exon to include XCE entry (.bam input only)"
         arg_type = Int
         default  = 1000
      "--aberrant-xri-one-side"
         help     = "Allow XRI calls where only one exon-intron spanning junction has support (.bam input only)"
         action   = :store_true
      "--no-aberrant"
         help     = "Don't allow aberrant splicing calls (.bam input only)"
         action   = :store_false

      "--sam"
         help     = "Should SAM format be sent to stdout? (.fastq input only)"
         action   = :store_true
      "--stranded"
         help     = "Is the data strand specific in fwd orientation? If so, increase speed with this flag (.fastq input only)"
         action   = :store_true
      "--pair-same-strand"
         help     = "Whippet by default tries to align fwd/rev pairs, if your data is fwd/fwd or rev/rev set this flag (.fastq input only)"
         action   = :store_true
      "--phred-33"
         help     = "Qual string is encoded in Phred+33 integers [default] (.fastq input only)"
         action   = :store_true
      "--phred-64"
         help     = "Qual string is encoded in Phred+64 integers"
         action   = :store_true
      "--circ"
         help     = "Allow back/circular splicing, this will allow output of `BS`-type lines (.fastq input only)"
         action   = :store_true
      "--force-gz"
         help     = "Regardless of suffix, consider read input as gzipped (.fastq input only)"
         action   = :store_true

      "--seed-len", "-L"
         help     = "Seed length (.fastq input only)"
         arg_type = Int
         default  = 18
      "--seed-try", "-M"
         help     = "Number of failed seeds to try before giving up (.fastq input only)"
         arg_type = Int
         default  = 3
      "--seed-tol", "-T"
         help     = "Number of seed hits to tolerate (.fastq input only)"
         arg_type = Int
         default  = 4
      "--seed-buf", "-B"
         help     = "Ignore this many bases from beginning and end of read for seed (.fastq input only)"
         arg_type = Int
         default  = 5
      "--seed-inc", "-I"
         help     = "Number of bases to increment seed each iteration (.fastq input only)"
         arg_type = Int
         default  = 18
      "--pair-range", "-P"
         help     = "Seeds for paired end reads must match within _ bases of one another (.fastq input only)"
         arg_type = Int
         default  = 5000
      "--mismatches", "-X"
         help     = "Allowable number of mismatches in alignment [counted as 1-10^(-phred/10)] (.fastq input only)"
         arg_type = Int
         default  = 3
      "--score-min", "-S"
         help     = "Minimum percent matching (matches - mismatches) / read_length (.fastq input only)"
         arg_type = Float64
         default  = 0.7

   end
   return parse_args(s)
end

function _resolve_batch_input(raw::AbstractString, manifest_dir::AbstractString)
   raw_s = String(raw)
   p = fixpath(raw_s)
   isfile(p) && return p
   if !isabspath(raw_s)
      alt = fixpath(joinpath(manifest_dir, raw_s))
      isfile(alt) && return alt
   end
   return p
end

function parse_batch_manifest(manifest_path::AbstractString)
   mpath = fixpath(String(manifest_path))
   isfile(mpath) || error("Batch manifest file not found: $mpath")
   mdir = dirname(abspath(mpath))
   entries = Tuple{String,String,Union{Nothing,String},String}[]
   delim = nothing
   saw_first = false

   open(mpath, "r") do io
      for (line_no, raw_line) in enumerate(eachline(io))
         line = strip(raw_line)
         if isempty(line) || startswith(line, '#')
            continue
         end
         if delim === nothing
            delim = occursin('\t', line) ? '\t' : ','
         end
         cols = [String(strip(c)) for c in split(line, delim)]
         if !saw_first
            saw_first = true
            if lowercase(cols[1]) in ("sample", "sample_id")
               continue
            end
         end
         if length(cols) < 3
            error("Manifest line $line_no needs at least 3 columns (sample,r1,out_prefix or sample,r1,r2,out_prefix): $line")
         end
         sname = cols[1]
         r1 = _resolve_batch_input(cols[2], mdir)
         if length(cols) >= 4
            r2_raw = cols[3]
            r2 = (isempty(r2_raw) || r2_raw == "-" || uppercase(r2_raw) == "NA") ? nothing : _resolve_batch_input(r2_raw, mdir)
            out_pref = fixpath(cols[4])
         else
            r2 = nothing
            out_pref = fixpath(cols[3])
         end
         push!(entries, (sname, r1, r2, out_pref))
      end
   end
   isempty(entries) && error("Batch manifest contains no sample rows: $mpath")
   return entries
end

function quantify_one_sample( lib, args, indexpath::String, inputfile::String, matefile::Union{Nothing,String}, outprefix::String )
   isbaminput = hasextension(inputfile, ".bam")

   if isbaminput
      ispaired = ispaired_bamfile(inputfile)
   else
      ispaired = matefile !== nothing
   end

   ContainerType = ispaired ? SGAlignPaired : SGAlignSingle

   if args["biascorrect"]
      CounterType   = JointBiasCounter
      BiasModelType = JointBiasMod
   else
      CounterType   = DefaultCounter
      BiasModelType = DefaultBiasMod
   end

   param = AlignParam( args, ispaired, kmer=lib.kmer )
   quant = GraphLibQuant{ContainerType,CounterType}( lib )
   multi = MultiMapping{ContainerType,CounterType}()
   mod   = BiasModelType()

   enc_offset = args["phred-64"] ? 64 : 33

   outdir = dirname(abspath(outprefix))
   isdir(outdir) || mkpath(outdir)

   println(stderr, "Processing reads from file...")
   if Threads.nthreads() > 1
      println(stderr, "Using $(Threads.nthreads()) threads for read alignment...")
   end
   mate_parser = nothing
   if isbaminput
      parser = ispaired ? make_collated_parser( inputfile ) : make_single_parser( inputfile )
      println(stderr, "BAM: " * inputfile)

      @timer concordant,nonconcordant,singletons,readlen = read_collated_bam( parser, lib, quant, mod )
      
      if ispaired
         println(stderr, "Finished loading $(concordant+nonconcordant+singletons) total alignments")
         println(stderr, "\tConcordant reads: $concordant")
         println(stderr, "\tNon-concordant reads: $nonconcordant")
         println(stderr, "\tSingle-mapping reads: $singletons")
      else
         println(stderr, "Finished loading $singletons valid alignments")
      end
      readlen  = Int(floor(readlen))
      mapped   = concordant + singletons
      totreads = concordant + nonconcordant + singletons

   else  # default fastq input
      parser = make_fqparser( inputfile, forcegzip=args["force-gz"] )

      if ispaired
         mate_parser = make_fqparser( matefile, forcegzip=args["force-gz"] )
         println(stderr, "FASTQ_1: " * inputfile)
         println(stderr, "FASTQ_2: " * matefile)
         @timer mapped,totreads,readlen = process_paired_reads!( parser, mate_parser, param, lib, quant, multi, mod,
                                                           sam=args["sam"], qualoffset=enc_offset )
         readlen = Int(floor(readlen))
         println(stderr, "Finished mapping $mapped paired-end reads of length $readlen each out of a total of $totreads mate-pairs...")
      else
         println(stderr, "FASTQ: " * inputfile)
         @timer mapped,totreads,readlen = process_reads!( parser, param, lib, quant, multi, mod,
                                                     sam=args["sam"], qualoffset=enc_offset )
         readlen = Int(floor(readlen))
         println(stderr, "Finished mapping $mapped single-end reads of length $readlen out of a total of $totreads reads...")
      end
   end
   close(parser)
   mate_parser !== nothing && close(mate_parser)

   # TPM_EM
   println(stderr, "Calculating expression values and MLE of equivalence classes with EM:")
   primer_normalize!( mod )
   primer_adjust!( quant, mod )
   primer_adjust!( multi, mod )
   println(stderr, "- $( length(multi.map) ) multi-gene mapping read equivalence classes...")
   build_equivalence_classes!( quant, lib, assign_long=true )
   println(stderr, "- $( length(quant.classes) ) multi-isoform equivalence classes...")
   calculate_tpm!( quant, readlen=readlen )
   @timer iter = gene_em!( quant, multi, sig=1, readlen=readlen, maxit=10000 )
   println(stderr, "Finished calculating transcripts per million (TpM) after $iter iterations of EM...")
   set_gene_tpm!( quant, lib )
   gc_normalize!( mod, lib, quant )
   gc_adjust!( quant, mod )
   gc_adjust!( multi, mod )

   output_tpm( outprefix, lib, quant )
   output_stats( outprefix * ".map.gz", lib, quant, param, indexpath, totreads, mapped, Int(round(total_multi(multi))), readlen, ver )
   output_junctions( outprefix * ".jnc.gz", lib, quant )

   println(stderr, "Assigning multi-mapping reads based on maximum likelihood estimate..")
   # Now assign multi to edges.
   @timer assign_ambig!( quant, lib, multi )

   effective_lengths!( lib, quant, 1, 0)
   println(stderr, "Calculating maximum likelihood estimate of events..." )
   @timer process_events( outprefix * ".psi.gz" , lib, quant, iscircok=args["circ"], readlen=readlen )
end

function main()

   args = parse_cmd()

   println(stderr, " $( round( (time_ns()-start)/1e9, digits=6 ) ) seconds." )

   batch_file = strip(args["batch"])
   has_batch  = !isempty(batch_file)
   has_input  = args["filename.fastq[.gz]--OR--pre-aligned.bam"] !== nothing

   if has_batch && has_input
      error("Specify either positional input FASTQ/BAM OR --batch <manifest.csv>, not both.")
   elseif !has_batch && !has_input
      error("Must provide an input FASTQ/BAM file or --batch <manifest.csv>.")
   end
   if has_batch && args["sam"]
      error("--sam cannot be combined with --batch (stdout streams would interleave across samples).")
   end

   indexpath = fixpath( args["index"] )
   indexname = hasextension(indexpath, ".jls") ? indexpath : indexpath * ".jls"
   println(stderr, "Loading splice graph index... $indexname")
   @timer lib = open(deserialize, indexname)
   checkversion( lib, ver, minversion )
   lib = FastGraphLib( lib )

   if has_batch
      entries = parse_batch_manifest(batch_file)
      nsamp = length(entries)
      println(stderr, "Batch mode: quantifying $nsamp sample(s) from $batch_file (index loaded once)...")
      for (si, (sname, r1, r2, out_pref)) in enumerate(entries)
         println(stderr, "\n=== [$si/$nsamp] Sample: $sname -> $out_pref ===")
         quantify_one_sample( lib, args, indexpath, r1, r2, out_pref )
         if si < nsamp
            GC.gc()
         end
      end
      println(stderr, "\nWhippet $ver batch done ($nsamp samples)." )
   else
      inputfile = fixpath( args["filename.fastq[.gz]--OR--pre-aligned.bam"] )
      matefile  = args["paired_mate.fastq[.gz]"] === nothing ? nothing : fixpath(args["paired_mate.fastq[.gz]"])
      outprefix = fixpath( args["out"] )
      quantify_one_sample( lib, args, indexpath, inputfile, matefile, outprefix )
      println(stderr, "Whippet $ver done." )
   end
end

@timer main()
