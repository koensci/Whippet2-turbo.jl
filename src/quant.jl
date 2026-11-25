
const SCALING_FACTOR = 1_000_000

# use specialty type for "read counts"
# with basic function: 'sum' that can be
# overridden for more complex count bias models
#const ReadCount = JointBiasCounter #DefaultCounter
const DEFAULT_DEF   = (v = zero(DefaultCounter); v.isadjusted = true; v)
const DEFAULT_JOINT = (v = zero(JointBiasCounter); v.isadjusted = true; v)

const DEF_READCOUNT = 0.0
const DEF_READVALUE = 1.0

default(::Type{DefaultCounter}) = DEFAULT_DEF
default(::Type{JointBiasCounter}) = DEFAULT_JOINT

# Use this abstraction to store SGAlignPaths before
# quantification or assignment
abstract type SGAlignContainer end

struct SGAlignSingle <: SGAlignContainer
   fwd::SGAlignPath
end

struct SGAlignPaired <: SGAlignContainer
   fwd::SGAlignPath
   rev::SGAlignPath
end

SGAlignPaired( single_path::SGAlignPath ) = SGAlignPaired( single_path, single_path )

# Store unique path's for gene/node but ignore
# SGAlignScore which is irrelevant at this stage
Base.isequal( a::SGAlignSingle, b::SGAlignSingle ) = a.fwd == b.fwd
Base.isequal( a::SGAlignPaired, b::SGAlignPaired ) = a.fwd == b.fwd && a.rev == b.rev

Base.:(==)( a::SGAlignSingle, b::SGAlignSingle ) = a.fwd == b.fwd
Base.:(==)( a::SGAlignPaired, b::SGAlignPaired ) = a.fwd == b.fwd && a.rev == b.rev

Base.hash( v::SGAlignSingle ) = hash( v.fwd )
Base.hash( v::SGAlignSingle, h::UInt64 ) = hash( v.fwd, h )
Base.hash( v::SGAlignPaired ) = v.fwd == v.rev ? 
                                 hash( v.fwd ) : 
                                 hash( v.fwd, hash( v.rev ) )
Base.hash( v::SGAlignPaired, h::UInt64 ) = v.fwd == v.rev ? 
                                 hash( v.fwd, h ) : 
                                 hash( v.fwd, hash( v.rev, h ) )

ispaired( cont::SGAlignSingle ) = false
ispaired( cont::SGAlignPaired ) = true

isexonic( node::SGAlignNode )  = true
isexonic( node::SGNodeIsExon ) = node.meta

isaberrant( node::NodeFloat ) = true
isaberrant( node::NodeInt )   = false

function unique_push!( arr::Vector{T}, el::K ) where {T,K}
   tel = convert(T, el)
   if !( tel in arr )
      push!( arr, tel )
   end
end

# check if one alignment path is a subset of another
function Base.in( first::SGAlignPath, last::SGAlignPath )
   for i in 1:(length(last)-length(first))+1
      offset = i - 1
      match = true
      for j in 1:length(first)
         if !(first[j] == last[offset+j])
            match = false
            break
         end
      end
      if match
         return true
      end
   end
   return false
end


# reimplementation of to edge code in events.jl
# we require not just the nodes to exist in the intset
# but also that they are adjacent in the intset
function Base.in( first, last, is::S ) where S <: MonotonicSet
   if first in is &&
      last in is
      for i in first+1:last-1
         if i in is
            return false
         end
      end
      return true
   end
   return false
end

Base.in( edge::I, is::BitSet ) where I <: IntervalValue = in( edge.first, edge.last, is )
Base.in( edge::IntervalValue{K,V}, is::SortedSet{K, O} ) where {K <: NodeFloat, V, O <: Base.Ordering} = in( edge.first, edge.last, is )


# This looks for an edge in a Vector of BitSets using ^
function Base.in( edge::I, viset::Vector{S} ) where {I <: IntervalValue, S <: MonotonicSet}
   for iset in viset
      if edge in iset
         return true
      end
   end
   false
end

# For events.jl:
function inall( edge::I, viset::Vector{S} ) where {I <: IntervalValue, S <: MonotonicSet}
   inone = false
   for iset in viset
      if edge.first >= first(iset) && edge.last <= last(iset)
         if !(edge in iset)
            return false
         end
         inone = true
      end
   end
   inone ? true : false
end

# each edge in the path must exist in the int set to be true
function Base.in( path::SGAlignPath, is::S ) where S <: MonotonicSet
   if length(path) == 1
      if !(path[1].node in is)
         return false
      end
   else
      for i in 1:length(path)-1
         if !in( path[i].node, path[i+1].node, is )
            return false
         end
      end
   end
   return true
end

Base.in( aln::SGAlignSingle, is::S ) where S <: MonotonicSet = in( aln.fwd, is )
Base.in( aln::SGAlignPaired, is::S ) where S <: MonotonicSet = in( aln.fwd, is ) && in( aln.rev, is )

# This is where we count reads for nodes/edges/circular-edges/effective_lengths
# bias is an adjusting multiplier for nodecounts to bring them to the same level
# as junction counts, which should always be at a lower level
mutable struct SpliceGraphQuant{C <: SGAlignContainer, R <: ReadCounter}
   node::Vector{R}
   edge::IntervalMap{NodeInt,R}
   aber::IntervalMap{NodeFloat,R}
   lnod::Vector{NodeFloat}
   rnod::Vector{NodeFloat}
   long::Dict{C,R}
   circ::Dict{Tuple{ExonInt,ExonInt},R}
   leng::Vector{Float64}
   bias::Float64

   function SpliceGraphQuant{C,R}( sg::AnySpliceGraph ) where {C <: SGAlignContainer, R <: ReadCounter}
      return new( R[zero(R) for x in 1:length(sg.nodelen)],
             IntervalMap{NodeInt,R}(),
             IntervalMap{NodeFloat,R}(),
             Vector{NodeFloat}(),
             Vector{NodeFloat}(),
             Dict{C,R}(),
             Dict{Tuple{ExonInt,ExonInt},R}(),
             zeros( length(sg.nodelen) ), 
             1.0 )
   end
end

abstract type EquivalenceClass end

const ClassInt  = Int32
const ClassType = Vector{ClassInt}

# type-stable isoform compatibility class
mutable struct IsoCompat <: EquivalenceClass
   class::ClassType
   prop::Vector{Float64}
   prop_sum::Float64
   count::Float64
   isdone::Bool

   function IsoCompat( arr::Vector{Int}, count::Float64 )
      return new( arr, ones( length(arr) ) / length(arr), 1.0, count, false )
   end
end

# Here we store whole graphome quantification
struct GraphLibQuant{C <: SGAlignContainer, R <: ReadCounter}
   tpm::Vector{Float64}                 # isoform tpm
   count::Vector{Float64}               # isoform counts
   length::Vector{Int64}                # isoform lengths
   geneidx::Vector{Int64}               # 0-based offset of isoforms for gene
   genetpm::Vector{Float64}             # sum of isoform tpms for each gene
   quant::Vector{SpliceGraphQuant{C,R}} # splice graph quant structs for gene
   classes::Vector{IsoCompat}           # isoform compatibility classes

   function GraphLibQuant{C,R}( lib::SeqLibrary ) where {C <: SGAlignContainer, R <: ReadCounter}
      ngenes = length(lib.graphs)
      isonum = Vector{Int}(undef, ngenes)
      isolen = 0
      @inbounds for i in 1:ngenes
         ni = length(lib.graphs[i].annoname)
         isonum[i] = ni
         isolen += ni
      end
      tpm     = zeros( Float64, isolen )
      count   = zeros( Float64, isolen )
      leng    = zeros( Int64, isolen )
      geneoff = Vector{Int64}(undef, ngenes)
      genetpm = zeros( Float64, ngenes )
      quant   = Vector{SpliceGraphQuant{C,R}}(undef, ngenes)
      classes = Vector{IsoCompat}()
      cumul_i = 0
      @inbounds for i in 1:ngenes
         geneoff[i] = cumul_i
         sg = lib.graphs[i]
         nlen = sg.nodelen
         apaths = sg.annopath
         for j in 1:length(apaths)
            s = Int64(0)
            for y in apaths[j]
               s += Int64(nlen[y])
            end
            leng[cumul_i + j] = s
         end
         cumul_i += isonum[i]
         quant[i] = SpliceGraphQuant{C,R}( sg )
      end
      new( tpm, count, leng, geneoff, genetpm, quant, classes )
   end
end

# Adjust each MultiCompat raw value and assign its adjusted to .count
function adjust!( sgquant::SpliceGraphQuant{C,R}, mod::B, func=adjust! ) where {C <: SGAlignContainer, R <: ReadCounter, B <: BiasModel}
   for node in sgquant.node
      func( node, mod )
   end
   for edg in sgquant.edge
      func( edg.value, mod )
   end
   for edg in sgquant.aber
      func( edg.value, mod )
   end
   for lng in values(sgquant.long)
      func( lng, mod )
   end
   for cir in values(sgquant.circ)
      func( cir, mod )
   end
   sgquant
end

primer_adjust!( sgquant::SpliceGraphQuant{C,R}, mod::B ) where {C <: SGAlignContainer, R <: ReadCounter, B <: BiasModel} = adjust!( sgquant, mod, primer_adjust! )
gc_adjust!( sgquant::SpliceGraphQuant{C,R}, mod::B ) where {C <: SGAlignContainer, R <: ReadCounter, B <: BiasModel} = adjust!( sgquant, mod, gc_adjust! )

# Adjust each sgquant storage value
function adjust!( gquant::GraphLibQuant{C,R}, mod::B, func=adjust! ) where {C <: SGAlignContainer, R <: ReadCounter, B <: BiasModel}
   for sgq in gquant.quant
      func( sgq, mod )
   end
   gquant
end

primer_adjust!( gquant::GraphLibQuant{C,R}, mod::B ) where {C <: SGAlignContainer, R <: ReadCounter, B <: BiasModel} = adjust!( gquant, mod, primer_adjust! )
gc_adjust!( gquant::GraphLibQuant{C,R}, mod::B ) where {C <: SGAlignContainer, R <: ReadCounter, B <: BiasModel} = adjust!( gquant, mod, gc_adjust! )

function gc_normalize!( mod::JointBiasMod, lib::SeqLibrary, quant::GraphLibQuant )
   @inbounds for i in 1:length(lib.graphs)
      sg  = lib.graphs[i]
      idx = quant.geneidx[i]
      for j in 1:length(sg.annobias)
         w = quant.tpm[j+idx]
         w > 0.0 && add_to_back!( mod.gcback, sg.annobias[j], w )
      end
   end
   foresum = sum(mod.gcfore)
   backsum = sum(mod.gcback)
   @inbounds for i in 1:length(mod.gcfore)
      if foresum > 0.0
         @fastmath mod.gcfore[i] /= foresum
      end
      if backsum > 0.0
         @fastmath mod.gcback[i] /= backsum
      end
   end
   mod
end

@inbounds function set_gene_tpm!( quant::GraphLibQuant, lib::SeqLibrary )
   for i in 1:length(lib.graphs)
      tpm = 0.0
      curidx = quant.geneidx[i]
      for j in 1:length( lib.graphs[i].annoname )
         tpm += quant.tpm[ curidx + j ]
      end
      quant.genetpm[i] = tpm
   end
end

function countedge!( sgquant::SpliceGraphQuant{C,R}, 
                     a::SGNodeMeta{N,A}, 
                     b::SGNodeMeta{M,A},
                     value ) where {N <: NodeNum, M <: NodeNum, A <: Any, C <: SGAlignContainer, R <: ReadCounter}

   lnode = a.node
   rnode = b.node
   if isaberrant(lnode) || isaberrant(rnode)
      interv = IntervalTrees.Interval{NodeFloat}( lnode, rnode )
      pushzero!( sgquant.aber, interv, value )
      isaberrant(lnode) && unique_push!( a.meta ? sgquant.lnod : sgquant.rnod, lnode )
      isaberrant(rnode) && unique_push!( b.meta ? sgquant.lnod : sgquant.rnod, rnode )
   elseif lnode < rnode # canonical splicing
      interv = IntervalTrees.Interval{NodeInt}( lnode, rnode )
      pushzero!( sgquant.edge, interv, value )
   elseif lnode >= rnode # back splicing
      pushzero!( sgquant.circ, (lnode,rnode), value )
   end
end

assign_path!( graphq::GraphLibQuant{SGAlignSingle}, 
              path::SGAlignSingle,
              value, 
              temp_iset::BitSet ) = assign_path!( graphq, path, value )

# This function re-count!'s SGAlignPaths which were not count!'ed originally
function assign_path!( graphq::GraphLibQuant{SGAlignSingle}, 
                       path::SGAlignSingle, 
                       value,
                       temp_iset::S ) where S <: MonotonicSet
   assign_path!( graphq, path.fwd, value )
end

function assign_path!( graphq::GraphLibQuant{SGAlignSingle}, 
                       path::SGAlignSingle, 
                       value)
   assign_path!( graphq, path.fwd, value )
end

@inbounds function assign_path!( graphq::GraphLibQuant{A}, 
                                 path::SGAlignPath, 
                                 value ) where A <: Union{SGAlignSingle, SGAlignPaired}
   init_gene = path[1].gene
   sgquant = graphq.quant[ init_gene ]

   if length(path) == 1
      if !isaberrant(path[1].node)
         push!( sgquant.node[ path[1].node ], value )
      end
   else
      for i in 1:length(path)-1
         countedge!( sgquant, path[i], path[i+1], value )
      end
   end
end

# This function re-counts paired SGAlignPaths that weren't counted because they were too long, disjointed, or both
@inbounds function assign_path!( graphq::GraphLibQuant{SGAlignPaired}, 
                                 path::SGAlignPaired, 
                                 value, 
                                 temp_iset::S=SortedSet{NodeFloat}() ) where S <: MonotonicSet
   if path.fwd == path.rev
      assign_path!( graphq, path.fwd, value )
      return
   end

   init_gene = path.fwd[1].gene
   sgquant = graphq.quant[ init_gene ]

   if length(path.fwd) == 1
      if !isaberrant(path.fwd[1].node)
         push!( sgquant.node[ path.fwd[1].node ], value )
         push!( temp_iset, path.fwd[1].node )
      end
   else
      for i in 1:length(path.fwd)-1
         lnode = path.fwd[i].node
         rnode = path.fwd[i+1].node
         push!(temp_iset, lnode)
         push!(temp_iset, rnode)

         countedge!( sgquant, path.fwd[i], path.fwd[i+1], value )
      end
   end

   if length(path.rev) == 1
      if !(path.rev[1].node in temp_iset) && !isaberrant(path.rev[1].node)
         push!( sgquant.node[ path.rev[1].node ], value )
      end
   else
      for i in 1:length(path.rev)-1
         lnode = path.rev[i].node
         rnode = path.rev[i+1].node
         (lnode in temp_iset && rnode in temp_iset) && continue

         countedge!( sgquant, path.rev[i], path.rev[i+1], value )
      end
   end
   empty!(temp_iset)
end

# build compatibility classes from SpliceGraphQuant node and edge "counts"
@inbounds function build_equivalence_classes!( graphq::GraphLibQuant, lib::SeqLibrary; assign_long=true )
   # initialize re-used data structures
   sset = SortedSet{NodeFloat}()
   iset = BitSet()
   temp = Dict{BitSet,Float64}()

   @inline function handle_iset_value!( idx, value )
      n_iso = length(iset)
      if n_iso == 1
         graphq.count[first(iset) + idx] += value
      elseif n_iso > 1
         ht_idx = Base.ht_keyindex(temp, iset)
         if ht_idx > 0
            @inbounds temp.vals[ht_idx] += value
         else
            temp[copy(iset)] = value
         end
      end
   end

   # go through splice graphs (to re-use temp data structures)
   for i in 1:length(lib.graphs)
      sg = lib.graphs[i]
      sgq = graphq.quant[i]
      apaths = sg.annopath
      n_anno = length(apaths)
      gidx = graphq.geneidx[i]

      # go through single node mapping reads
      nodes = sgq.node
      for n in 1:length(nodes)
         val = get(nodes[n])
         val > 0.0 || continue
         for p in 1:n_anno
            if n in apaths[p]
               push!(iset, p)
            end
         end
         handle_iset_value!( gidx, val )
         empty!(iset)
      end
      # go through edges
      for edg in sgq.edge
         val = get(edg.value)
         val > 0.0 || continue
         for p in 1:n_anno
            if edg in apaths[p]
               push!(iset, p)
            end
         end
         handle_iset_value!( gidx, val )
         empty!(iset)
      end
      # go through multi-edge crossing paths
      for (align, cnt) in sgq.long
         val = get(cnt)
         if val > 0.0
            for p in 1:n_anno
               if align in apaths[p]
                  push!(iset, p)
               end
            end
            handle_iset_value!( gidx, val )
            empty!(iset)
         end
         if assign_long
            empty!(sset)
            assign_path!( graphq, align, val, sset )
         end
      end
      # create compatibility classes for isoforms
      if !isempty(temp)
         for (bset, val) in temp
            if val > 0.0
               arr = Vector{Int}(undef, length(bset))
               k = 1
               for p in bset
                  arr[k] = p + gidx
                  k += 1
               end
               push!( graphq.classes, IsoCompat(arr, val) )
            end
         end
         empty!(temp)
      end
   end
end


# This function takes a vector of alignments and returns the compatibility class
# and a boolean describing whether or not we can assign the read through EM
@inbounds function equivalence_class( graphq::GraphLibQuant,
                                      lib::SeqLibrary,
                                      cont::Vector{C},
                                      temp_iset::BitSet=BitSet() ) where C <: SGAlignContainer
   matches_all = true
   for path in cont
      g = path.fwd[1].gene
      has_match = false
      apaths = lib.graphs[g].annopath
      gidx = graphq.geneidx[g]
      for p in 1:length(apaths)
         if path in apaths[p]
            push!(temp_iset, p + gidx)
            has_match = true
         end
      end
      if !has_match
         matches_all = false
         empty!(temp_iset)
         return ClassType(),false
      end
   end
   retval = collect(ClassInt, temp_iset)
   empty!(temp_iset)
   retval,matches_all
end


## MultiAlignment containers & Equivalence classes

mutable struct MultiCompat{R} <: EquivalenceClass
   class::ClassType
   prop::Vector{Float64}
   prop_sum::Float64
   count::Float64
   isdone::Bool
   postassign::Bool
   raw::R

   function MultiCompat{R}( graphq::GraphLibQuant{C,R}, lib::SeqLibrary, cont::Vector{C}, temp_iset::BitSet=BitSet() ) where {C <: SGAlignContainer, R <: ReadCounter}
      class,matches = equivalence_class( graphq, lib, cont, temp_iset )
      new( class, ones(length(class)) / length(class), 1.0, 1.0, !matches, !matches, zero(R) )
   end
end

# This hash structure stores multi-mapping
# equivalence classes
struct MultiMapping{C <: SGAlignContainer, R <: ReadCounter}
   map::Dict{Vector{C},MultiCompat{R}}
   iset::BitSet

   function MultiMapping{C,R}() where {C <: SGAlignContainer, R <: ReadCounter}
      new( Dict{Vector{C},MultiCompat{R}}(), BitSet() )
   end
end

function total_multi( mm::MultiMapping )
   total = 0.0
   for v in values(mm.map)
      total += v.count
   end
   total
end

function Base.push!( ambig::MultiMapping{SGAlignSingle,R}, 
                     alns::Vector{SGAlignment}, 
                     value,
                     graphq::GraphLibQuant, 
                     lib::SeqLibrary ) where R <: ReadCounter
   cont = Vector{SGAlignSingle}(undef, length(alns))
   @inbounds for i in 1:length(alns)
      cont[i] = SGAlignSingle( alns[i].path )
   end
   ht_idx = Base.ht_keyindex( ambig.map, cont )
   if ht_idx > 0
      @inbounds push!( ambig.map.vals[ht_idx].raw, value )
   else
      @inbounds for i in 1:length(cont)
         cont[i] = SGAlignSingle( copy(alns[i].path) )
      end
      mc = MultiCompat{R}( graphq, lib, cont, ambig.iset )
      push!( mc.raw, value )
      ambig.map[cont] = mc
   end
end

function Base.push!( ambig::MultiMapping{SGAlignPaired,R}, 
                     fwd::Vector{SGAlignment}, 
                     rev::Vector{SGAlignment}, 
                     value,
                     graphq::GraphLibQuant, 
                     lib::SeqLibrary ) where R <: ReadCounter
   cont = Vector{SGAlignPaired}( undef, length(fwd) )
   @inbounds for i in 1:length(fwd)
      cont[i] = SGAlignPaired( fwd[i].path, rev[i].path )
   end
   ht_idx = Base.ht_keyindex( ambig.map, cont )
   if ht_idx > 0
      @inbounds push!( ambig.map.vals[ht_idx].raw, value )
   else
      @inbounds for i in 1:length(cont)
         cont[i] = SGAlignPaired( copy(fwd[i].path), copy(rev[i].path) )
      end
      mc = MultiCompat{R}( graphq, lib, cont, ambig.iset )
      push!( mc.raw, value )
      ambig.map[cont] = mc
   end
end

# Adjust each MultiCompat raw value and assign its adjusted to .count
function adjust!( multi::MultiMapping{C,R}, 
                  mod::B, 
                  func=adjust! ) where {C <: SGAlignContainer, 
                                        R <: ReadCounter, 
                                        B <: BiasModel}
   for mc in values(multi.map)
      func( mc.raw, mod )
      mc.count = get(mc.raw)
   end
   multi
end

primer_adjust!( multi::MultiMapping{C,R}, mod::B ) where {C <: SGAlignContainer, R <: ReadCounter, B <: BiasModel} = adjust!( multi, mod, primer_adjust! )
gc_adjust!( multi::MultiMapping{C,R}, mod::B ) where {C <: SGAlignContainer, R <: ReadCounter, B <: BiasModel} = adjust!( multi, mod, gc_adjust! )

# for a given alignment container, set the temp_iset in multi to the compatibility class
@inbounds function set_equivalence_class!( multi::MultiMapping{C,R}, 
                                           graphq::GraphLibQuant,
                                           lib::SeqLibrary, 
                                           aln::C ) where {C <: SGAlignContainer, R <: ReadCounter}
   empty!(multi.iset)
   g = aln.fwd[1].gene
   apaths = lib.graphs[g].annopath
   gidx = graphq.geneidx[g]
   for p in 1:length(apaths)
      if aln in apaths[p]
         push!(multi.iset, p + gidx)
      end
   end
end

function get_class_value( iset::BitSet, compat::MultiCompat )
   value = 0.0
   @inbounds for i in 1:length(compat.class)
      if compat.class[i] in iset
         value += compat.prop[i]
      end
   end
   value
end

function assign_ambig!( graphq::GraphLibQuant{C,R}, 
                        lib::SeqLibrary, 
                        ambig::MultiMapping{C,R} ) where {C <: SGAlignContainer, R <: ReadCounter}
   vals = Float64[]
   @inbounds for (vp,mc) in ambig.map
      n_vp = length(vp)
      resize!(vals, n_vp)

      if mc.postassign  # not all paths match annotated transcripts, use gene level values to assign
         total_tpm = 0.0
         for i in 1:n_vp
            v = graphq.genetpm[ vp[i].fwd[1].gene ]
            vals[i] = v
            total_tpm += v
         end
         if total_tpm > 0.0
            for i in 1:n_vp
               @fastmath vals[i] /= total_tpm
            end
         else
            inv_n = 1.0 / n_vp
            for i in 1:n_vp
               vals[i] = inv_n
            end
         end
      else # assign based on transcripts
         total_v = 0.0
         for i in 1:n_vp
            set_equivalence_class!( ambig, graphq, lib, vp[i] )
            v = get_class_value( ambig.iset, mc )
            vals[i] = v
            total_v += v
         end
         if total_v > 0.0
            for i in 1:n_vp
               @fastmath vals[i] /= total_v
            end
         else
            inv_n = 1.0 / n_vp
            for i in 1:n_vp
               vals[i] = inv_n
            end
         end
      end

      # now assign alignment paths proportional to their weight
      for i in 1:n_vp
         assign_path!( graphq, vp[i], vals[i] * mc.count, ambig.iset )
      end
   end
end

# Count a single alignment in quant
function count!( graphq::GraphLibQuant{C,R}, 
                 align::SGAlignment, 
                 value ) where {C <: SGAlignContainer, R <: ReadCounter}

   align.isvalid == true || return
   count!( graphq, align.path, value )
end

function count!( graphq::GraphLibQuant{C,R}, 
                 path::SGAlignPath, 
                 value ) where {C <: SGAlignContainer, R <: ReadCounter}

   init_gene = path[1].gene
   sgquant   = graphq.quant[ init_gene ]

   # single node support
   if length(path) == 1
      # access node -> [ SGNode( gene, *node* ) ]
      !isaberrant(path[1].node) && push!( sgquant.node[ path[1].node ], value )
   else # multi-node mapping read
      # if any node in path maps to other gene, ignore read
      for n in 1:length(path)
         if path[n].gene != init_gene
            return
         end
      end

      # add a single edge support
      if length(path) == 2
         countedge!( sgquant, path[1], path[2], value )
      else # or we have long path, add long count
         curkey = C(path)
         ht_idx = Base.ht_keyindex( sgquant.long, curkey )
         if ht_idx > 0
            @inbounds push!( sgquant.long.vals[ht_idx], value )
         else
            entry = zero(R)
            push!( entry, value )
            sgquant.long[C(copy(path))] = entry
         end
      end
   end
end

## Extension of count! for paired end counting
function count!( graphq::GraphLibQuant{C,R}, 
                 fwd::SGAlignment, 
                 rev::SGAlignment, 
                 value ) where {C <: SGAlignContainer, R <: ReadCounter}
   
   (fwd.isvalid == true && rev.isvalid == true) || return
   count!( graphq, fwd.path, rev.path, value )
end

function count!( graphq::GraphLibQuant{C,R}, 
                 fwd::SGAlignPath, 
                 rev::SGAlignPath, 
                 value ) where {C <: SGAlignContainer, R <: ReadCounter}

   init_gene = fwd[1].gene
   rev_gene  = rev[1].gene
   (init_gene == rev_gene) || return
   sgquant = graphq.quant[ init_gene ]
   singlebool = false

   if length(fwd) <= length(rev)
      if fwd in rev
         singlebool = true
         single = rev
      end
   else
      if rev in fwd
         singlebool = true
         single = fwd
      end
   end

   # single node support
   if singlebool
       count!( graphq, single, value )
   else
       # or we have two paths, add long count
       curkey = SGAlignPaired(fwd, rev)
       ht_idx = Base.ht_keyindex( sgquant.long, curkey )
       if ht_idx > 0
          @inbounds push!( sgquant.long.vals[ht_idx], value )
       else
          entry = zero(R)
          push!( entry, value )
          sgquant.long[SGAlignPaired(copy(fwd), copy(rev))] = entry
       end
   end
end

@inline function calculate_tpm!( quant::GraphLibQuant, 
                                 counts::Vector{Float64}=quant.count; 
                                 readlen::Int64=50, 
                                 sig::Int64=1 )
   tpm = quant.tpm
   leng = quant.length
   rpk_sum = 0.0
   @inbounds for i in 1:length(counts)
      c = counts[i]
      if c != 0.0
         v = @fastmath c / max( Float64(leng[i] - readlen), 1.0 )
         tpm[i] = v
         rpk_sum += v
      else
         tpm[i] = 0.0
      end
   end
   rpk_sum = max( rpk_sum, 1.0 )
   scale = @fastmath SCALING_FACTOR / rpk_sum
   if sig > 0
      @inbounds for i in 1:length(tpm)
         v = tpm[i]
         if v != 0.0
            @fastmath tpm[i] = round( v * scale, digits=sig )
         end
      end
   else
      @inbounds for i in 1:length(tpm)
         v = tpm[i]
         if v != 0.0
            @fastmath tpm[i] = v * scale
         end
      end
   end
end

calculate_bias!( sgquant::SpliceGraphQuant ) = 1.0
global_bias( graphq::GraphLibQuant ) = 1.0,0.0

# Every exon-exon junction has an effective mappable space of readlength - K*2.
# Since we only count nodes that don't map over edges, we can give
# give each junction an effective_length of one, and now the space
# inside of a node for which a read could map without overlapping a junction
# is put into units of junction derived effective_length
# kadj is minimum of (readlength - minalignlen) and k-1
@inline function eff_length( node, sg::AnySpliceGraph, eff_len::Int, kadj::Int )
   len = sg.nodelen[node] + (istxstart( sg.edgetype[node] ) ? 0 : kadj) +
                            (istxstop( sg.edgetype[node+1] ) ? 0 : kadj)
   @fastmath len / eff_len
end

function eff_lengths!( sg::AnySpliceGraph, sgquant::SpliceGraphQuant, eff_len::Int, kadj::Int )
   @inbounds for i in 1:length( sg.nodelen )
      sgquant.leng[i] = sg.nodelen[i] #eff_length( i, sg, eff_len, kadj )
   end
end

function effective_lengths!( lib::SeqLibrary, graphq::GraphLibQuant, eff_len::Int, kadj::Int )
   @inbounds for i in 1:length( lib.graphs )
      eff_lengths!( lib.graphs[i], graphq.quant[i], eff_len, kadj )
   end
end

function unsafe_copy!( dest::Vector{T}, src::Vector{T}; indx_shift=0 ) where T <: Number
   @inbounds for i in 1:length(src)
      dest[i+indx_shift] = src[i]
   end
end

function copy_isdone!( dest::Vector{T}, src::Vector{T}, denominator::Float64=1.0, sig::Int=4 ) where T <: AbstractFloat
   if denominator == 0.0
      @inbounds for i in 1:length(dest)
         dest[i] = zero(T)
      end
      return true
   end
   isdifferent = false
   @inbounds for i in 1:length(dest)
      s = src[i]
      val = s == 0.0 ? zero(T) : round( s / denominator, digits=sig )
      if val != dest[i]
         isdifferent = true
      end
      dest[i] = isnan(val) ? zero(T) : val
   end
   !isdifferent
end

# This function performs expectation maximization
# and sets the quant.tpm array as the proposed expression set
# at each iteration.   It also requires that calculate_tpm! be
# run initially once prior to gene_em!() call.

function gene_em!( quant::GraphLibQuant, 
                   ambig::MultiMapping{C,R};
                   it::Int64=1, 
                   maxit::Int64=1000,
                   sig::Int64=4, 
                   readlen::Int64=50 ) where {C <: SGAlignContainer, R <: ReadCounter}

   n_iso      = length(quant.count)
   count_temp = copy(quant.count)
   tpm        = quant.tpm
   qcount     = quant.count
   qleng      = quant.length

   # Collect active (unconverged) equivalence classes into dense vectors and track active isoforms
   is_active_iso = falses(n_iso)
   @inbounds for i in 1:n_iso
      if qcount[i] != 0.0 || tpm[i] != 0.0
         is_active_iso[i] = true
      end
   end

   max_class_len = 500
   active_ambig = MultiCompat{R}[]
   sizehint!(active_ambig, length(ambig.map))
   for eq in values(ambig.map)
      eq.isdone && continue
      push!(active_ambig, eq)
      cl = eq.class
      length(cl) > max_class_len && (max_class_len = length(cl))
      @inbounds for idx in cl
         is_active_iso[idx] = true
      end
   end

   active_classes = IsoCompat[]
   sizehint!(active_classes, length(quant.classes))
   for eq in quant.classes
      eq.isdone && continue
      push!(active_classes, eq)
      cl = eq.class
      length(cl) > max_class_len && (max_class_len = length(cl))
      @inbounds for idx in cl
         is_active_iso[idx] = true
      end
   end

   active_iso = Int[]
   sizehint!(active_iso, count(is_active_iso))
   @inbounds for i in 1:n_iso
      is_active_iso[i] && push!(active_iso, i)
   end

   # Precompute inverse effective lengths for active isoforms
   n_act_iso = length(active_iso)
   inv_eff_len = Vector{Float64}(undef, n_act_iso)
   @inbounds for k in 1:n_act_iso
      i = active_iso[k]
      inv_eff_len[k] = 1.0 / max( Float64(qleng[i] - readlen), 1.0 )
   end

   prop_temp = zeros(Float64, max_class_len)

   @inline function maximize_and_assign!( eq::E, maximize::Bool ) where E <: EquivalenceClass
      cl = eq.class
      pr = eq.prop
      cnt = eq.count
      n_cl = length(cl)
      if maximize
         psum = 0.0
         @inbounds for c in 1:n_cl
            cur_tpm = tpm[cl[c]]
            prop_temp[c] = cur_tpm
            psum += cur_tpm
         end
         eq.prop_sum = psum
         eq.isdone = copy_isdone!( pr, prop_temp, psum )
      end
      if eq.isdone
         @inbounds for c in 1:n_cl
            idx = cl[c]
            inc = pr[c] * cnt
            qcount[idx] += inc
            count_temp[idx] += inc
         end
      else
         @inbounds for c in 1:n_cl
            count_temp[cl[c]] += pr[c] * cnt
         end
      end
      return eq.isdone
   end

   changed = true
   while changed && it < maxit
      @inbounds for k in 1:n_act_iso
         i = active_iso[k]
         count_temp[i] = qcount[i]
      end

      do_max = (it > 1)

      # Process active multi-mapping classes and compact out newly converged ones
      w = 0
      @inbounds for j in 1:length(active_ambig)
         eq = active_ambig[j]
         if !maximize_and_assign!( eq, do_max )
            w += 1
            active_ambig[w] = eq
         end
      end
      resize!(active_ambig, w)

      # Process active isoform classes and compact out newly converged ones
      w = 0
      @inbounds for j in 1:length(active_classes)
         eq = active_classes[j]
         if !maximize_and_assign!( eq, do_max )
            w += 1
            active_classes[w] = eq
         end
      end
      resize!(active_classes, w)

      # Update TPM only for active isoforms and check convergence inline
      rpk_sum = 0.0
      @inbounds for k in 1:n_act_iso
         i = active_iso[k]
         v = @fastmath count_temp[i] * inv_eff_len[k]
         count_temp[i] = v
         rpk_sum += v
      end
      rpk_sum = max( rpk_sum, 1.0 )
      scale = @fastmath SCALING_FACTOR / rpk_sum

      changed = false
      if sig > 0
         @inbounds for k in 1:n_act_iso
            i = active_iso[k]
            v = count_temp[i]
            new_tpm = v == 0.0 ? 0.0 : round( @fastmath(v * scale), digits=sig )
            if new_tpm != tpm[i]
               changed = true
               tpm[i] = new_tpm
            end
         end
      else
         @inbounds for k in 1:n_act_iso
            i = active_iso[k]
            v = count_temp[i]
            new_tpm = @fastmath v * scale
            if new_tpm != tpm[i]
               changed = true
               tpm[i] = new_tpm
            end
         end
      end

      it += 1
   end
   it
end


function output_tpm( file, lib::SeqLibrary, gquant::GraphLibQuant )

   geneio = open( file * ".gene.tpm.gz", "w" )
   genestream = ZlibDeflateOutputStream( geneio )
   isoio = open( file * ".isoform.tpm.gz", "w" )
   isostream = ZlibDeflateOutputStream( isoio )
   output_tpm_header( genestream, "Gene" )
   output_tpm_header( isostream, "Isoform" )
   @inbounds for i in 1:length(lib.names)
      idx = gquant.geneidx[i]
      genetpm = 0.0
      genecnt = 0.0
      annonames = lib.graphs[i].annoname
      for j in 1:length(annonames)
         tval = gquant.tpm[j+idx]
         cval = gquant.count[j+idx]
         genetpm += tval
         genecnt += cval
         write( isostream, annonames[j], UInt8('\t') )
         write_float( isostream, tval )
         write( isostream, UInt8('\t') )
         write_float( isostream, cval )
         write( isostream, "\t\n" )
      end
      write( genestream, lib.names[i], UInt8('\t') )
      write_float( genestream, genetpm )
      write( genestream, UInt8('\t') )
      write_float( genestream, genecnt )
      write( genestream, UInt8('\n') )
   end
   close( isostream )
   close( isoio )
   close( genestream )
   close( geneio )
end
