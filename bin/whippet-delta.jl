#!/usr/bin/env julia
# Tim Sterne-Weiler 2016
# Upgraded with Beta-Binomial GLM, Cox-Reid Empirical Bayes dispersion shrinkage,
# arbitrary covariate formulas, LRT p-values, and Benjamini-Hochberg FDR.

using Pkg

const dir = abspath( splitdir(@__FILE__)[1] )
const ver = chomp(readline(open(dir * "/VERSION")))

start = time_ns()
println( stderr, "Whippet $ver loading... " )

Pkg.activate(dir * "/..")

using Whippet
using Random
using ArgParse
using Glob

function parse_cmd()
   s = ArgParseSettings(
      description = "Differential splicing analysis using a Beta-Binomial GLM with Cox-Reid " *
                    "Empirical Bayes dispersion shrinkage, covariate support, and BH FDR " *
                    "(or original Whippet Beta Monte Carlo via --legacy)."
   )
   @add_arg_table! s begin
    "--a", "-a"
      help     = "Replicates for Set A -- glob pattern matching *.psi* or comma-delimited list of filenames (2-group mode)."
      arg_type = String
    "--b", "-b"
      help     = "Replicates for Set B -- same rules as (-a) (2-group mode)."
      arg_type = String
    "--design"
      help     = "Path to a CSV/TSV sample metadata table for GLM covariate modeling. Contains a file path column (--file-col, default 'file') plus predictor/covariate columns."
      arg_type = String
    "--file-col"
      help     = "Column name in --design containing paths to .psi/.psi.gz files."
      arg_type = String
      default  = "file"
    "--formula", "-f"
      help     = "DESeq2-style model formula for --design mode, e.g. '~ batch + donor + condition' or '~ batch * condition'."
      arg_type = String
    "--test"
      help     = "Term or contrast to test in --design mode: e.g. 'condition', 'condition:treated_vs_ctrl', or a continuous column like 'age'. Defaults to the last term in --formula."
      arg_type = String
      default  = ""
    "--out", "-o"
      help     = "Core file name prefix to send .diff.gz output to!"
      arg_type = String
      default  = fixpath( "./output" )
    "--directory", "-d"
      help     = "Directory to search for file patterns or relative file paths in -a/-b or --design."
      arg_type = String
      default  = "."
    "--min-reads", "-r"
      help     = "Minimum number of reads for an event in a sample to be included."
      arg_type = Int64
      default  = 5
    "--min-samples", "-s"
      help     = "Minimum number of valid samples (per group in 2-group/contrast mode, or total in general GLM mode)."
      arg_type = Int64
      default  = 1
    "--no-shrink"
      help     = "Disable Empirical Bayes dispersion shrinkage (use per-event Cox-Reid dispersion MLE)."
      action   = :store_true
    "--no-stratify-eb"
      help     = "Disable separate Empirical Bayes dispersion priors for spliced (CE/AA/AD/RI/AF/AL) vs. tandem UTR (TS/TE) events."
      action   = :store_true
    "--fixed-dispersion"
      help     = "Use a fixed Beta-Binomial overdispersion phi in (0, 1) instead of estimating from data (disabled if <= 0)."
      arg_type = Float64
      default  = -1.0
    "--legacy"
      help     = "Use legacy Whippet Beta Monte Carlo differential testing (only valid with -a and -b)."
      action   = :store_true
    "--pseudo-adjust"
      help     = "[Legacy mode] When a single replicate is given, pseudo-adjust a \"second\" for MLE fitting by this value."
      arg_type = Float64
      default  = 0.01
    "--max-variance"
      help     = "[Legacy mode] Maximum allowed variance across psi values to fit a Beta -- if greater, use null Beta(0.1,0.1)."
      arg_type = Float64
      default  = 0.15
    "--use-depth"
      help     = "[Legacy mode] Sample each replicate's PSI according to read-depth of the event."
      action   = :store_true
    "--emperical-size", "-e"
      help     = "[Legacy mode] Empirical distribution size to sample from."
      arg_type = Int64
      default  = 1000
    "--seed", "-g"
      help     = "Seed the RNG (Int) for reproducible results on successive runs."
      arg_type = Int64
      default  = 123456
   end
   return parse_args(s)
end

function retrievefilelist( pattern::String, dir::String )
   list = Vector{String}()
   if something(findfirst(isequal(','), pattern), 0) > 0
      tmp = split( pattern, ',', keepempty=false )
      for f in tmp
         fs = strip(String(f))
         push!( list, isabspath(fs) ? fs : joinpath(dir, fs) )
      end
   else
      tmp = glob( "*" * pattern * "*.psi*", dir )
      for file in tmp
         push!( list, string(file) )
      end
   end
   list
end

function main()
   args  = parse_cmd()
   println(stderr, " $( round( (time_ns()-start)/1e9, digits=6 ) ) seconds" )
   Random.seed!( args["seed"] )
   dir   = fixpath( args["directory"] )
   outfile = args["out"] * ".diff.gz"

   has_design = args["design"] !== nothing
   has_ab     = args["a"] !== nothing || args["b"] !== nothing

   if has_design && has_ab
      error("Specify either --design (with --formula) OR -a/-b, not both.")
   end
   if !has_design && !has_ab
      error("Must specify either -a and -b (2-group comparison) OR --design and --formula (covariate GLM).")
   end

   if args["legacy"]
      has_design && error("--legacy mode only supports -a and -b (covariates require the default Beta-Binomial GLM engine).")
      (args["a"] === nothing || args["b"] === nothing) && error("Both -a and -b are required.")
      lista = retrievefilelist( args["a"], dir )
      listb = retrievefilelist( args["b"], dir )
      println(stderr, "Sample A: $(join(map(basename, lista), ','))")
      println(stderr, "Sample B: $(join(map(basename, listb), ','))")
      if length(lista) <= 0 || length(listb) <= 0
         error("Unable to match files! length(a) == $(length(lista)), length(b) == $(length(listb))!")
      end
      astreams = open_streams( lista )
      bstreams = open_streams( listb )
      println(stderr, "Now processing files using legacy Beta Monte Carlo method...")
      @timer process_psi_files( outfile, astreams, bstreams,
                                min_samp=args["min-samples"],
                                min_reads=args["min-reads"],
                                amt=0.0,
                                size=args["emperical-size"],
                                point_est=!args["use-depth"],
                                pseudo_adj=args["pseudo-adjust"],
                                max_variance=args["max-variance"] )
   else
      if has_design
         args["formula"] === nothing && error("--formula is required when using --design (e.g. --formula '~ batch + condition').")
         design = build_design_from_table( fixpath(args["design"]), args["formula"];
                                           test_spec=args["test"], dir=dir,
                                           file_col=args["file-col"] )
         println(stderr, "Design table: $(args["design"]) ($(length(design.files)) samples)")
         println(stderr, "Formula:      $(args["formula"])")
         println(stderr, "Testing:      $(design.test_label) (columns: $(join(design.col_names[design.test_cols], ", ")))")
      else
         (args["a"] === nothing || args["b"] === nothing) && error("Both -a and -b are required.")
         lista = retrievefilelist( args["a"], dir )
         listb = retrievefilelist( args["b"], dir )
         println(stderr, "Sample A: $(join(map(basename, lista), ','))")
         println(stderr, "Sample B: $(join(map(basename, listb), ','))")
         if length(lista) <= 0 || length(listb) <= 0
            error("Unable to match files! length(a) == $(length(lista)), length(b) == $(length(listb))!")
         end
         design = build_two_group_design( lista, listb )
      end

      println(stderr, "Running Beta-Binomial GLM differential splicing analysis...")
      @timer process_psi_files_glm( outfile, design;
                                    min_samp=args["min-samples"],
                                    min_reads=args["min-reads"],
                                    shrink_dispersion=!args["no-shrink"],
                                    stratify_eb=!args["no-stratify-eb"],
                                    fixed_phi=args["fixed-dispersion"] )
   end

   println(stderr, "Whippet $ver done." )
end

@timer main()
