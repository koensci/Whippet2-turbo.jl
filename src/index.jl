# index.jl - tim.sterne.weiler@utoronto.ca, 1/28/16


abstract type SeqLibrary end

struct GraphLib <: SeqLibrary
   offset::Vector{CoordInt}
   names::Vector{GeneName}
   info::Vector{GeneInfo}
   lengths::Vector{Float64}
   graphs::Vector{SpliceGraph}
   coords::IntervalCollection{SGNodeIsExon}
   edges::Edges
   index::FMIndex
   sorted::Bool
   kmer::Int64
   hasindex::Bool             
   version::VersionNumber
end

function _has_ambig_dna4( seq::LongSequence{DNAAlphabet{4}} )
   len = length(seq)
   len == 0 && return false
   if first(seq.part) == 1
      nfull = len >>> 4
      data = seq.data
      @inbounds for w_idx in 1:nfull
         w = data[w_idx]
         a = w & 0x1111111111111111
         b = (w >> 1) & 0x1111111111111111
         c = (w >> 2) & 0x1111111111111111
         d = (w >> 3) & 0x1111111111111111
         (a + b + c + d) != 0x1111111111111111 && return true
      end
      @inbounds for i in ((nfull << 4) + 1):len
         isambiguous(seq[i]) && return true
      end
      return false
   else
      @inbounds for i in 1:len
         isambiguous(seq[i]) && return true
      end
      return false
   end
end

# Concretely-typed zero-copy in-memory wrappers for hot paths.
# Preserves 100% serialization compatibility with existing .jls GraphLib files.
struct FastSpliceGraph{K}
   nodeoffset::Vector{CoordInt}
   nodecoord::Vector{CoordInt}
   nodelen::Vector{CoordInt}
   edgetype::Vector{EdgeType}
   edgeleft::Vector{SGKmer{K}}
   edgeright::Vector{SGKmer{K}}
   annopath::Vector{BitSet}
   annoname::Vector{String}
   annobias::Vector{ExpectedGC}
   seq::LongSequence{DNAAlphabet{4}}
   has_ambig::Bool
end

@inline _sg_has_ambig( sg::FastSpliceGraph ) = sg.has_ambig
@inline _sg_has_ambig( sg::SpliceGraph ) = true

struct FastGraphLib{K, F <: FMIndex} <: SeqLibrary
   offset::Vector{CoordInt}
   names::Vector{GeneName}
   info::Vector{GeneInfo}
   lengths::Vector{Float64}
   graphs::Vector{FastSpliceGraph{K}}
   coords::IntervalCollection{SGNodeIsExon}
   edges::Edges{K}
   index::F
   sorted::Bool
   kmer::Int64
   hasindex::Bool
   version::VersionNumber
end

const AnySpliceGraph = Union{SpliceGraph, FastSpliceGraph}

function FastGraphLib( lib::GraphLib )
   K = Int(lib.kmer)
   edges = lib.edges::Edges{K}
   n = length(lib.graphs)
   graphs = Vector{FastSpliceGraph{K}}(undef, n)
   @inbounds for i in 1:n
      g = lib.graphs[i]::SpliceGraph{K}
      seq = g.seq isa LongSequence{DNAAlphabet{4}} ? g.seq : LongSequence{DNAAlphabet{4}}(g.seq)
      graphs[i] = FastSpliceGraph{K}(g.nodeoffset, g.nodecoord, g.nodelen, g.edgetype,
                                     g.edgeleft, g.edgeright, g.annopath, g.annoname,
                                     g.annobias, seq, _has_ambig_dna4(seq))
   end
   return FastGraphLib{K, typeof(lib.index)}(
      lib.offset, lib.names, lib.info, lib.lengths, graphs,
      lib.coords, edges, lib.index, lib.sorted, lib.kmer, lib.hasindex, lib.version
   )
end

FastGraphLib( lib::FastGraphLib ) = lib


function checkversion( lib::GraphLib, currentver, minversion::VersionNumber )
   try
      if lib.version < minversion
         error("Index is from Whippet $(lib.version)\nWhippet $currentver requires a $minversion index or greater!")
      end
   catch e
      error("Index is from Whippet 1.6 or older!\nWhippet $currentver requires a $minversion index or greater!")
   end
   true
end

function build_chrom_dict( ref::RefSet )
   ret = Dict{SeqName,Vector{GeneName}}() # refgenomeseq->geneid[]
   for g in keys(ref)
      chrom = ref[g].info.name
      if !haskey(ret, chrom)
         ret[chrom] = Vector{GeneName}()
      end
      push!(ret[chrom], g)
   end
   ret
end

function twobit_enc(seq)
   len = length(seq)
   ret = IntVector{2,UInt8}(len)
   @inbounds for i in 1:len
      b = seq[i]
      val = isambiguous(b) ? 0x00 : convert(UInt8, trailing_zeros(b))
      ret[i] = val
   end
   ret
end

function twobit_enc_graphs( graphs::Vector{SpliceGraph}, total_len::Int )
   ret = IntVector{2,UInt8}(total_len)
   pos = 1
   @inbounds for g in graphs
      seq = g.seq::LongSequence{DNAAlphabet{4}}
      len = length(seq)
      for i in 1:len
         b = seq[i]
         ret[pos] = isambiguous(b) ? 0x00 : convert(UInt8, trailing_zeros(b))
         pos += 1
      end
   end
   return ret
end


function load_fasta( fhIter; verbose=false )
   seq = SGSequence(mutable=false)
   offset = Int[]
   names = SeqName[]
   for r in fhIter
      immutable!(r.seq)
      push!(names, r.name)
      push!(offset, length(seq)+1)
      println( stderr, r )
      @time seq *= r.seq
   end
   seq, offset, names
end

function single_genome_index!( fhIter; verbose=false )
   seq,offset,names = load_fasta( fhIter )
   @time fm = FMIndex(twobit_enc(seq), 4, r=4, program=:SuffixArrays, mmap=true)
   println( stderr, "Finished building Index..." )

   SeqLib(seq, offset, names, fm, false)
end


function trans_index!( fhIter, ref::RefSet; kmer=9, version=v"1.7" )
   seqdic  = build_chrom_dict( ref )
   ng_est  = length(ref)
   xoffset = sizehint!(Vector{CoordInt}(), ng_est)
   xgenes  = sizehint!(Vector{GeneName}(), ng_est)
   xinfo   = sizehint!(Vector{GeneInfo}(), ng_est)
   xlength = sizehint!(Vector{Float64}(), ng_est)
   xgraph  = sizehint!(Vector{SpliceGraph}(), ng_est)
   # set up splice graphs (parallel per chromosome, deterministic order, zero O(N^2) xcript concat)
   runoffset = 0
   num = 0
   for r in fhIter
      name = String(r.data[r.identifier])
      haskey(seqdic, name) || continue
      sg = LongSequence{DNAAlphabet{4}}( r.data[r.sequence] )

      genes = seqdic[name]
      ng = length(genes)
      println(stderr, "Building Splice Graphs for $( name ).." )
      @timer begin
         chr_graphs = Vector{SpliceGraph}(undef, ng)
         if Threads.nthreads() > 1 && ng > 1
            Threads.@threads for gi in 1:ng
               @inbounds chr_graphs[gi] = SpliceGraph( ref[genes[gi]], sg, kmer )
            end
         else
            @inbounds for gi in 1:ng
               chr_graphs[gi] = SpliceGraph( ref[genes[gi]], sg, kmer )
            end
         end
         @inbounds for gi in 1:ng
            g = genes[gi]
            curgraph = chr_graphs[gi]
            rg = ref[g]
            push!(xgraph, curgraph)
            push!(xgenes, g)
            push!(xinfo, rg.info)
            push!(xlength, rg.length)
            push!(xoffset, CoordInt(runoffset))
            runoffset += length(curgraph.seq)
            num += 1
         end
      end
   end
   if num > 1
      println( stderr, "Building full sg-index from $num genes..." )
      @time fm = FMIndex(twobit_enc_graphs(xgraph, runoffset), 4, r=1, program=:SuffixArrays, mmap=true)
   else
      error("ERROR: No genes in the GTF file matched chromosome names in the FASTA file!")
   end

   GC.gc()

   println( stderr, "Building edges.." )
   @time edges = build_edges( xgraph, kmer ) # TODO make variable kmer
   @time ic = interval_index( xgraph, xinfo, min_intron_size=MININTRONSIZE )

   GraphLib( xoffset, xgenes, xinfo, xlength, xgraph, ic, edges, fm, true, kmer, true, version )
end

function fasta_to_index( filename::String, ref::RefSet; kmer=9 )
   if isgzipped( filename )
      println(stderr, "Decompressing and Indexing $filename...")
      to_open = open( filename ) |> x->ZlibInflateInputStream(x, reset_on_end=true, bufsize=262144)
   else
      println(stderr, "Indexing $filename...")
      to_open = open( filename ) |> x->BufferedInputStream(x, 262144)
   end
   # iterate through fasta entries
   index = @time trans_index!(FASTA.Reader( to_open ), ref, kmer=kmer)
   index
end

function add_nodes_to_collection!( ic::IntervalCollection{SGNodeIsExon}, 
                                   info::GeneInfo, 
                                   sg::SpliceGraph, 
                                   gene::Int;
                                   min_intron_size=25 )
   for i in 1:length(sg.nodecoord)
      sg.nodelen[i] <= 0 && continue
      left     = sg.nodecoord[i]
      right    = sg.nodecoord[i]+sg.nodelen[i]-1
      exon = GenomicFeatures.Interval(info.name, 
                                      left, 
                                      right, 
                                      info.strand ? '+' : '-', 
                                      SGNodeIsExon(gene, NodeInt(i) :: NodeNum, true))
      push!( ic, exon )

      i < length(sg.nodecoord) || break

      if info.strand
         nextleft = sg.nodecoord[i+1]

         if sg.edgetype[i+1] == EDGETYPE_LR &&
            right + min_intron_size < nextleft

            intron = GenomicFeatures.Interval(info.name, 
                                             right+1, 
                                             nextleft-1, 
                                             '+', 
                                             SGNodeIsExon(gene, NodeFloat(i+0.5) :: NodeNum, false))
            push!( ic, intron )
         end
      else
         left  = sg.nodecoord[i+1]+sg.nodelen[i+1]
         right = sg.nodecoord[i]-1

         if sg.edgetype[i+1] == EDGETYPE_LR &&
            left + min_intron_size < right

            intron = GenomicFeatures.Interval(info.name, 
                                             left, 
                                             right, 
                                             '-', 
                                             SGNodeIsExon(gene, NodeFloat(i+0.5) :: NodeNum, false))
            push!( ic, intron )
         end
      end

   end
   ic
end

function interval_index( graphs::Vector{SpliceGraph},
                         infos::Vector{GeneInfo};
                         min_intron_size=25 )
   
   ic = IntervalCollection{SGNodeIsExon}()

   for i in 1:length(graphs)
      add_nodes_to_collection!(ic, infos[i], graphs[i], i, min_intron_size=min_intron_size)
   end
   ic
end
