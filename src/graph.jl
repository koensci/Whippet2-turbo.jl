# requires:

primitive type EdgeType 8 end

const EDGETYPE_TO_UINT8 = fill( 0x07, (4,4) )
      EDGETYPE_TO_UINT8[4,2] = 0x00 # 'SL' = 0x00; Tx Start
      EDGETYPE_TO_UINT8[4,3] = 0x01 # 'SR' = 0x01; Tandem Last Exon
      EDGETYPE_TO_UINT8[2,4] = 0x02 # 'LS' = 0x02; Tandem First Exon
      EDGETYPE_TO_UINT8[3,4] = 0x03 # 'RS' = 0x03; Tx End
      EDGETYPE_TO_UINT8[2,3] = 0x04 # 'LR' = 0x04; Intron 5'->3' SS
      EDGETYPE_TO_UINT8[2,2] = 0x05 # 'LL' = 0x05; Alt- 5'SS
      EDGETYPE_TO_UINT8[3,3] = 0x06 # 'RR' = 0x06; Alt- 3'SS

const INDEX_TO_EDGETYPE = transpose(reshape([[0x03,0x01,0x03,0x03];
                                             [0x06  for _ in 1:4 ];
                                             [0x05,0x04,0x05,0x02];
                                             [0x00  for _ in 1:4 ] ], (4,4)))

const INDEX_TO_EDGETYPE_NODE = transpose(reshape([[0x03  for _ in 1:4 ];
                                                  [0x06  for _ in 1:4 ];
                                                  [0x05  for _ in 1:4 ];
                                                  [0x00  for _ in 1:4 ] ], (4,4)))


function Base.convert( ::Type{EdgeType}, one::UInt8, two::UInt8 )
   @assert( 5 <= one <= 7 && 5 <= two <= 7 )
   EDGETYPE_TO_UINT8[one-3,two-3]
end

Base.convert( ::Type{EdgeType}, one::DNA, two::DNA ) =
              Base.convert( EdgeType, convert(UInt8, trailing_zeros(one)), convert(UInt8, trailing_zeros(two)) )
Base.convert( ::Type{EdgeType}, edge::UInt8 ) = reinterpret(EdgeType, edge)
Base.convert( ::Type{UInt8}, edge::EdgeType ) = reinterpret(UInt8, edge)
Base.convert( ::Type{I}, edge::EdgeType) where {I <: Integer} = Base.convert(I, Base.convert(UInt8, edge))

EdgeType(edge::UInt8) = reinterpret(EdgeType, edge)
UInt8(edge::EdgeType) = reinterpret(UInt8, edge)

const EDGETYPE_SL = convert(EdgeType, 0b000)
const EDGETYPE_SR = convert(EdgeType, 0b001)
const EDGETYPE_LS = convert(EdgeType, 0b010)
const EDGETYPE_RS = convert(EdgeType, 0b011)

const EDGETYPE_LR = convert(EdgeType, 0b100)
const EDGETYPE_LL = convert(EdgeType, 0b101)
const EDGETYPE_RR = convert(EdgeType, 0b110)

istxstart(   edge::EdgeType ) = edge == EDGETYPE_SL || edge == EDGETYPE_LS ? true : false
istxstop(    edge::EdgeType ) = edge == EDGETYPE_SR || edge == EDGETYPE_RS ? true : false

isexonstart( edge::EdgeType ) = edge in (EDGETYPE_LR, EDGETYPE_RR, EDGETYPE_SR, EDGETYPE_LS, EDGETYPE_SL) ? true : false
isexonstop(  edge::EdgeType ) = edge in (EDGETYPE_LR, EDGETYPE_LL, EDGETYPE_SR, EDGETYPE_LS, EDGETYPE_RS) ? true : false
ishardedge(  edge::EdgeType ) = edge in (EDGETYPE_LR, EDGETYPE_LS, EDGETYPE_SR) ? true : false

isacceptor(  edge::EdgeType ) = edge in (EDGETYPE_LR, EDGETYPE_RR, EDGETYPE_SR) ? true : false
isdonor(     edge::EdgeType ) = edge in (EDGETYPE_LR, EDGETYPE_LL, EDGETYPE_LS) ? true : false


# Function takes coordinate types for node boundaries
# and returns an EdgeType
function get_edgetype( minidx::Int, secidx::Int, isnode::Bool, strand::Bool )
   if isnode
      ret = INDEX_TO_EDGETYPE_NODE[minidx,secidx]
   else
      ret = INDEX_TO_EDGETYPE[minidx,secidx]
   end
   ret = EdgeType(ret)
   if !strand
      ret = invert_edgetype( ret )
   end
   ret
end

# Invert EdgeTypes using xor
# 'LR' does not need to be inverted
function invert_edgetype( edge::EdgeType )
   if edge == EdgeType(0x04) # 'LR'
      return edge
   else
      return convert(EdgeType, convert(UInt8, edge) ⊻ 0b011)
   end
end

struct GCBiasParams
   length::Int64           # length offset
   width::Float64          # gc bin size
end

const ExpectedGC = Vector{Float16}

@inline _is_gc_dna( b::DNA ) = (b == DNA_G) | (b == DNA_C) | (b == DNA_S)

function ExpectedGC( seq::BioSequence{A}; gc_param = GCBiasParams(50, 0.05) ) where A <: BioSequences.Alphabet
   nbins = Int(div(1.0, gc_param.width) + 1)
   bins = zeros(Float16, nbins)
   wlen = Int(gc_param.length)
   slen = length(seq)
   slen < wlen && return bins

   counts = zeros(Int, nbins)
   inv_wlen = 1.0 / Float64(wlen)
   width = gc_param.width

   gc_cnt = 0
   @inbounds for i in 1:wlen
      gc_cnt += ifelse(_is_gc_dna(seq[i]), 1, 0)
   end
   @inbounds counts[Int(div(gc_cnt * inv_wlen, width)) + 1] += 1

   @inbounds for i in (wlen + 1):slen
      gc_cnt += ifelse(_is_gc_dna(seq[i]), 1, 0) - ifelse(_is_gc_dna(seq[i - wlen]), 1, 0)
      counts[Int(div(gc_cnt * inv_wlen, width)) + 1] += 1
   end

   total = slen - wlen + 1
   inv_tot = 1.0 / Float64(total)
   @inbounds for i in 1:nbins
      bins[i] = Float16(counts[i] * inv_tot)
   end
   return bins
end

# This holds a representation of the splice graph
# which is a directed multigraph
struct SpliceGraph{K}
   nodeoffset::Vector{CoordInt} # SG offset
   nodecoord::Vector{CoordInt}  # Genome offset
   nodelen::Vector{CoordInt}
   edgetype::Vector{EdgeType}
   edgeleft::Vector{SGKmer{K}}
   edgeright::Vector{SGKmer{K}}
   annopath::Vector{BitSet}
   annoname::Vector{String}
   annobias::Vector{ExpectedGC}
   seq::SGSequence
end
# All positive strand oriented sequences--->
# Node array: txStart| 1 |   2   | 3 |    4    |5| 6 |txEnd
# Edge array:        1   2       3   4         5 6   7
# Node coord:  chr   100 200     300 400      500 600 700

#= empty constructor
SpliceGraph(k::Int) = SpliceGraph( Vector{CoordInt}(), Vector{CoordInt}(),
                                   Vector{CoordInt}(), Vector{EdgeType}(),
                                   Vector{SGKmer{k}}(),Vector{SGKmer{k}}(), dna"" )
=#

# Main constructor
# Build splice graph here.
function SpliceGraph( gene::RefGene, genome::SGSequence, k::Int )
   # splice graph variables
   nodeoffset = Vector{CoordInt}()
   nodecoord  = Vector{CoordInt}()
   nodelen    = Vector{CoordInt}()
   edgetype   = Vector{EdgeType}()

   strand = gene.info.strand  # Bool

   # initialize iterators
   alen = length(gene.acc)
   dlen = length(gene.don)
   slen = length(gene.txst)
   plen = length(gene.txen)

   i1 = 1; i2 = 1; i3 = 1; i4 = 1

   # return a tuple containing the min coordinate's idx,val
   # each of 4 categories: 1=txen, 2=acc, 3=don, 4=txst
   @inline function getmin_ind_val( i1::Int, i2::Int, i3::Int, i4::Int )
      def = typemax(CoordInt)
      v1 = i1 <= plen ? (@inbounds gene.txen[i1]) : def
      v2 = i2 <= alen ? (@inbounds gene.acc[i2])  : def
      v3 = i3 <= dlen ? (@inbounds gene.don[i3])  : def
      v4 = i4 <= slen ? (@inbounds gene.txst[i4]) : def
      best_i = 1
      best_v = v1
      if v2 < best_v; best_i = 2; best_v = v2; end
      if v3 < best_v; best_i = 3; best_v = v3; end
      if v4 < best_v; best_i = 4; best_v = v4; end
      return best_i, best_v
   end

   while (i1 <= plen || i2 <= alen || i3 <= dlen || i4 <= slen)
      # iterate through donors, acceptors, txst, txen left to right;
      # '-' strand is reversed in O(N) after the loop.
      minidx, minval = getmin_ind_val(i1, i2, i3, i4)
      if minidx == 1; i1 += 1
      elseif minidx == 2; i2 += 1
      elseif minidx == 3; i3 += 1
      else; i4 += 1
      end
      secidx, secval = getmin_ind_val(i1, i2, i3, i4)

      # last coordinate in the gene:
      if secval == typemax(CoordInt)
         termedge = strand ? EdgeType(0x03) : EdgeType(0x00)
         push!(edgetype, termedge)
         break
      end

      # now should we make a node?
      if issubinterval( gene.exons, IntervalTrees.Interval{CoordInt}(minval,secval) ) ||
         (minidx == 2 && minval in gene.novelacc) ||
         (secidx == 3 && secval in gene.noveldon)
         leftadj  = (minidx == 1 || minidx == 3) && minval != secval ? 1 : 0
         rightadj = (secidx == 2 || secidx == 4) && minval != secval ? 1 : 0
         nodesize = Int(secval - minval) - leftadj - rightadj + 1
         edge     = get_edgetype( minidx, secidx, true, strand ) # determine EdgeType
         pushval  = minval + leftadj
      else # don't make a node, this is a sequence gap, make edge and inc+=2
         if secidx == 1; i1 += 1
         elseif secidx == 2; i2 += 1
         elseif secidx == 3; i3 += 1
         else; i4 += 1
         end
         thridx, thrval = getmin_ind_val(i1, i2, i3, i4)
         rightadj = (thridx == 2 || thridx == 4) && secval != thrval ? 1 : 0
         nodesize = Int(thrval - secval) - rightadj + 1
         edge     = get_edgetype( minidx, secidx, false, strand )
         pushval  = secval
      end

      push!(nodecoord, CoordInt(pushval))
      push!(nodelen,   CoordInt(max(0, nodesize)))
      push!(edgetype,  edge)

   end # end while

   if !strand
      reverse!(nodecoord)
      reverse!(nodelen)
      reverse!(edgetype)
   end

   # calculate node offsets and pre-allocate full gene sequence in O(N)
   nnodes = length(nodelen)
   sizehint!(nodeoffset, nnodes)
   curoffset = CoordInt(1)
   totlen = 0
   @inbounds for i in 1:nnodes
      push!(nodeoffset, curoffset)
      nl = Int(nodelen[i])
      curoffset += CoordInt(nl)
      totlen += nl
   end

   seq = LongSequence{DNAAlphabet{4}}(totlen)
   if strand
      @inbounds for i in 1:nnodes
         nl = Int(nodelen[i])
         nl > 0 || continue
         copyto!(seq, Int(nodeoffset[i]), genome, Int(nodecoord[i]), nl)
      end
   else
      @inbounds for i in 1:nnodes
         nl = Int(nodelen[i])
         nl > 0 || continue
         dst_off = Int(nodeoffset[i])
         src_end = Int(nodecoord[i]) + nl - 1
         for p in 0:(nl - 1)
            seq[dst_off + p] = BioSequences.complement(genome[src_end - p])
         end
      end
   end

   eleft  = Vector{SGKmer{k}}(undef, length(edgetype))
   eright = Vector{SGKmer{k}}(undef, length(edgetype))

   paths = build_paths_edges( nodecoord, nodelen, gene )
   names = map( x->x.info.name, gene.reftx )
   expgc = ExpectedGC[ExpectedGC(path_to_seq(p, nodeoffset, nodelen, seq)) for p in paths]

   return SpliceGraph{k}( nodeoffset, nodecoord, nodelen, edgetype, eleft, eright, paths, names, expgc, seq )
end

# re-orient - strand by using pushfirst! instead of push!
function stranded_push!( collection, value, strand::Bool )
   if strand
      push!( collection, value )
   else # '-' strand
      pushfirst!( collection, value )
   end
end

# This function looks specifically for an intersection
# between one interval, and the interval tree, such that
# some interval in the tree fully contains the sub interval
#  itree contains: (low,             high)
#  subint must be:     (low+x, high-y)
#  where x and y are >= 0
# Returns: bool
function issubinterval( itree::T, subint::I ) where {T <: IntervalTree,
                                                     I <: AbstractInterval}
   for i in intersect(itree, subint)
      if i.first <= subint.first && i.last >= subint.last
         return true
      end
   end
   false
end

# Note: Base.get(::Tuple, ::Integer, default) is provided natively by Julia Base.

# Functions for iso.jl annotated edges feature
# Take a RefTx and produce an BitSet through the nodes of a SpliceGraph
# returns: BitSet
function build_annotated_path( nodecoord::Vector{CoordInt},
                               nodelen::Vector{CoordInt},
                               tx::RefTx, 
                               strand::Bool )
   path = BitSet()
   nnodes = length(nodecoord)
   step = strand ? 1 : -1
   @inbounds for i in 1:length(tx.acc)
      rng = searchsorted( nodecoord, tx.acc[i], rev=!strand )
      ind = last(rng)
      push!( path, ind )
      cur = ind + step
      don_i = tx.don[i]
      while 1 <= cur <= nnodes && nodecoord[cur] <= don_i
         push!( path, cur )
         cur += step
      end
   end
   path
end

function build_paths_edges( nodecoord::Vector{CoordInt},
                            nodelen::Vector{CoordInt},
                            gene::RefGene )
   paths = Vector{BitSet}(undef, length(gene.reftx))
   @inbounds for (i, tx) in enumerate(gene.reftx)
      paths[i] = build_annotated_path( nodecoord, nodelen, tx, gene.info.strand )
   end
   paths
end

function path_to_seq( path::BitSet, 
                      nodeoffset::Vector{CoordInt},
                      nodelen::Vector{CoordInt}, 
                      seq::SGSequence )
   totlen = 0
   @inbounds for x in path
      totlen += Int(nodelen[x])
   end
   out = LongSequence{DNAAlphabet{4}}(totlen)
   pos = 1
   @inbounds for x in path
      nl = Int(nodelen[x])
      if nl > 0
         copyto!(out, pos, seq, Int(nodeoffset[x]), nl)
         pos += nl
      end
   end
   return out
end
