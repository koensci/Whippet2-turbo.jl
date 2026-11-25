
struct AlignParam
   mismatches::Int      # Allowable mismatches
   kmer_size::Int       # Minimum number of matches to extend past an edge
   seed_try::Int        # Starting number of seeds
   seed_tolerate::Int   # Allow at most _ hits for a valid seed
   seed_min_qual::Int   # Absolute lowest quality score to allow for worst base in seed
   seed_length::Int     # Seed Size
   seed_buffer::Int     # Ignore first and last _ bases
   seed_inc::Int        # Incrementation for subsequent seed searches
   seed_pair_range::Int # Seeds for paired-end reads must match within _ nt of each other
   score_range::Float64 # Identity scoring range to be considered repetitive
   score_min::Float64   # Minimum percent identity score for a valid alignment
   is_stranded::Bool    # Is input data + strand only?
   is_paired::Bool      # Paired end data?
   is_pair_rc::Bool     # Is mate 2 the reverse complement of mate 1
   is_trans_ok::Bool    # Do we allow edges from one gene to another
   is_circ_ok::Bool     # Do we allow back edges
end

AlignParam( ispaired = false ) = ispaired ? AlignParam( 3, 9, 2, 4, 4, 18, 5, 18, 2500, 0.05, 0.7,
                                                        false, true, true, false, true ) :
                                            AlignParam( 3, 9, 2, 4, 4, 18, 5, 18, 1000, 0.05, 0.7,
                                                        false, false, true, false, true )

# load from command line args
@inline function AlignParam( args::Dict{S,Any}, ispaired=false; kmer=9 ) where S <: AbstractString
   aln = AlignParam( args["mismatches"], kmer,
                           args["seed-try"],
                           args["seed-tol"], 
                           4,
                           args["seed-len"],
                           args["seed-buf"],
                           args["seed-inc"],
                           args["pair-range"], 
                           0.05,
                           args["score-min"],
                           args["stranded"], 
                           ispaired,
                           !args["pair-same-strand"], 
                           false,
                           args["circ"] )
   aln
end

abstract type UngappedAlignment end

const MatchesInt    = UInt8
const MismatchInt   = UInt8
const MismatchFloat = Float32

struct SGAlignScore
   matches::MatchesInt
   mismatches::MismatchInt
   mistolerance::MismatchFloat
end

@inline Base.zero(t::Type{SGAlignScore}) = SGAlignScore(0x0000,zero(MismatchInt),zero(MismatchFloat))
@inline Base.one(t::Type{SGAlignScore}) = SGAlignScore(0x0001,one(MismatchInt),one(MismatchFloat))
@inline match(t::Type{SGAlignScore}, val::I) where I <: Integer = SGAlignScore(convert(MatchesInt, val),zero(MismatchInt),zero(MismatchFloat))

const SGAlignNode = SGNodeMeta{NodeInt, SGAlignScore}

Base.hash( sgn::SGNodeMeta{J, A} ) where {J <: NodeNum, A <: Any} = hash( sgn.gene, hash(sgn.node) )
Base.hash( sgn::SGNodeMeta{J, A}, h::UInt64 ) where {J <: NodeNum, A <: Any} = hash( sgn.gene, hash(sgn.node, h) )
Base.isequal( a::SGNodeMeta{J, A}, b::SGNodeMeta{K, B} ) where {J <: NodeNum, A <: Any, K <: NodeNum, B <: Any} = a.gene == b.gene && 
                                                                                                                  a.node == b.node
Base.:(==)( a::SGNodeMeta{J, A}, b::SGNodeMeta{K, B} ) where {J <: NodeNum, A <: Any, K <: NodeNum, B <: Any} = a.gene == b.gene && 
                                                                                                                a.node == b.node

const SGAlignPath = Vector{SGNodeMeta{N, A}} where {N <: NodeNum, A <: Any}

@inline SGAlignNode(sgn::SGNode) = SGAlignNode(sgn.gene, sgn.node, zero(SGAlignScore))
@inline Base.push!( v::Vector{SGNodeMeta{N, A}}, sgn::SGNode ) where {N <: NodeNum, A <: Any} = push!( v, SGAlignNode(sgn) )

const ReadLengthInt = UInt8

mutable struct SGAlignment <: UngappedAlignment
   offset::UInt32
   offsetread::ReadLengthInt
   path::Vector{SGAlignNode}
   strand::Bool
   isvalid::Bool
end

SGAlignment() = SGAlignment(zero(UInt32), zero(ReadLengthInt), SGAlignNode[], true, false)

const SGAlignVec = Nullable{Vector{SGAlignment}}

const DEF_ALIGN = SGAlignment(zero(UInt32), zero(ReadLengthInt), SGAlignNode[], true, false)

@inbounds function matches( v::Vector{SGAlignNode} )
   matches = zero(MatchesInt)
   for i in 1:length(v)
      matches += v[i].meta.matches
   end
   matches
end
@inline @inbounds function matches_range( v::Vector{SGAlignNode}, a::Int, b::Int )
   m = zero(MatchesInt)
   for i in a:b
      m += v[i].meta.matches
   end
   m
end
@inbounds function mistolerance( v::Vector{SGAlignNode} )
   mistol = zero(MismatchFloat)
   for i in 1:length(v)
      mistol += v[i].meta.mistolerance
   end
   mistol
end
@inbounds function score( v::Vector{SGAlignNode} )
   score = zero(MismatchFloat)
   for i in 1:length(v)
      score += v[i].meta.matches - v[i].meta.mistolerance
   end
   score
end
@inline matches( aln::A ) where {A <: UngappedAlignment} = matches( aln.path )
@inline mistolerance( aln::A ) where {A <: UngappedAlignment} = mistolerance( aln.path )
@inline score( aln::A ) where {A <: UngappedAlignment} = score( aln.path )
@inline identity( align::A, readlen::I ) where {A <: UngappedAlignment, I <: Integer} = convert(Float32, @fastmath score( align ) / readlen)
@inline isvalid( align::A ) where {A <: UngappedAlignment} = align.isvalid && length(align.path) >= 1 && matches(align) > mistolerance(align) ? true : false

Base.:>( a::SGAlignment, b::SGAlignment ) = >( score(a), score(b) )
Base.:<( a::SGAlignment, b::SGAlignment ) = <( score(a), score(b) )
Base.:>=( a::SGAlignment, b::SGAlignment ) = >=( score(a), score(b) )
Base.:<=( a::SGAlignment, b::SGAlignment ) = <=( score(a), score(b) )
Base.isless( a::SGAlignment, b::SGAlignment ) = a < b

@inline copy_align( a::SGAlignment ) = SGAlignment(a.offset, a.offsetread, copy(a.path), a.strand, a.isvalid)

# Precomputed lookup table for phred -> probability of accurate base (1 - 10^(-phred/10))
const PHRED_TO_PROB_TABLE = Float32[convert(Float32, @fastmath 1 - 10^(-convert(Int8, reinterpret(Int8, UInt8(p))) / 10)) for p in 0:255]

# add prob of being accurate base to mismatch, rather than integer.
@inline phred_to_prob( phred::UInt8 ) = @inbounds PHRED_TO_PROB_TABLE[Int(phred) + 1]

@inline eq_dna2_dna4( rseq::LongSequence{DNAAlphabet{2}}, ridx::Int, sgseq::LongSequence{DNAAlphabet{4}}, sgidx::Int ) =
   (0x01 << raw_getindex(rseq, ridx)) == raw_getindex(sgseq, sgidx)
@inline eq_dna2_dna4( rseq, ridx::Int, sgseq, sgidx::Int ) = @inbounds rseq[ridx] == sgseq[sgidx]

@inline function Base.empty!( align::SGAlignment )
   align.offset = zero(UInt32)
   align.offsetread = zero(ReadLengthInt)
   empty!( align.path )
   align.strand = true
   align.isvalid = false
end

@inline function remove_invalid!( vec::Vector{SGAlignment} )
   i = 1
   while i <= length(vec)
      @inbounds if !isvalid( vec[i] )
         deleteat!( vec, i )
      else
         i += 1
      end
   end
end

@inline function seed_locate( p::AlignParam, 
                              index::FMIndex, 
                              read::FASTQRecord; 
                              offset_left::Bool=true, 
                              both_strands::Bool=true )
   def_sa = 2:1
   rseq = read.sequence
   rqual = read.quality
   readlen = length(rseq)
   seed_len = p.seed_length
   seed_tol = p.seed_tolerate
   seed_min_q = p.seed_min_qual
   seed_buf = p.seed_buffer
   ctry = 1
   if offset_left
      increment = p.seed_inc
      increment_sm = 1
      curpos = seed_buf
   else
      increment = -p.seed_inc
      increment_sm = -1
      curpos = readlen - seed_len - seed_buf
   end
   maxright = readlen - seed_len - seed_buf
   @inbounds while ctry <= p.seed_try && seed_buf <= curpos <= maxright
      qend = curpos + seed_len - 1
      low_qual = false
      for qidx in curpos:qend
         if rqual[qidx] <= seed_min_q
            low_qual = true
            break
         end
      end
      if low_qual
         curpos += increment_sm
         continue
      end
      sa = sa_range_sub( rseq, curpos, qend, index )
      cnt = length(sa)
      ctry += 1
      if cnt == 0 || cnt > seed_tol
         if both_strands
            rev_sa = sa_range_sub_rc( rseq, curpos, qend, index )
            rev_cnt = length(rev_sa)
            if 0 < rev_cnt <= seed_tol
               rev_pos = readlen - qend + 1
               return rev_sa, rev_pos, false
            end
         end
      else
         return sa, curpos, true
      end
      curpos += increment
   end
   def_sa, curpos, true
end

@inline function splice_by_score!( arr::Vector{A}, threshold, buffer ) where A <: UngappedAlignment
   i = 1
   while i <= length( arr )
      @inbounds if threshold - score( arr[i] ) > buffer
         splice!( arr, i )
         i -= 1
      end
      i += 1
   end
end

@inline function splice_by_identity!( arr::Vector{A}, threshold, buffer, readlen ) where A <: UngappedAlignment
   i = 1
   while i <= length( arr )
      @inbounds if threshold - identity( arr[i], readlen ) > buffer
         splice!( arr, i )
         i -= 1
      end
      i += 1
   end
end

@inline function reset_align!( align::SGAlignment, offset::UInt32, offsetread::ReadLengthInt, ispos::Bool )
   align.offset = offset
   align.offsetread = offsetread
   empty!( align.path )
   align.strand = ispos
   align.isvalid = false
   align
end

@inline function _ungapped_align!( p::AlignParam, 
                                   lib::SeqLibrary, 
                                   read::FASTQRecord, 
                                   indx::Int, 
                                   readloc::Int,
                                   align::SGAlignment;
                                   ispos::Bool=true, 
                                   geneind::GeneInt=convert(GeneInt,searchsortedlast( lib.offset, indx )) )

   sgidx       = indx - Int(@inbounds lib.offset[geneind])
   sg          = @inbounds lib.graphs[geneind]
   offset_node = convert(NodeInt, searchsortedlast(sg.nodeoffset, convert(CoordInt, sgidx + 1)))

   reset_align!( align, UInt32(sgidx + 1), convert(ReadLengthInt, readloc + 1), ispos )

   # If the splice graph has no ambiguous bases and the entire exact FM-index seed
   # stays inside offset_node, skip re-checking the (seed_length - 1) forward seed bases.
   node_end  = Int(@inbounds sg.nodeoffset[offset_node]) + Int(@inbounds sg.nodelen[offset_node]) - 1
   seed_skip = (!_sg_has_ambig(sg) && (sgidx + p.seed_length - 1 <= node_end) && (readloc + p.seed_length - 1 <= length(read.sequence))) ? (p.seed_length - 1) : 0

   align = ungapped_fwd_extend( p, lib, geneind, sgidx + 1 + seed_skip,
                                read, readloc + 1 + seed_skip, ispos=ispos,
                                align=align, nodeidx=offset_node,
                                init_mat=convert(MatchesInt, seed_skip) )

   align = ungapped_rev_extend( p, lib, geneind, sgidx,
                                read, readloc, ispos=ispos, align=align,
                                nodeidx=offset_node )
   align
end

@inline function _ungapped_align( p::AlignParam, 
                                  lib::SeqLibrary, 
                                  read::FASTQRecord, 
                                  indx::Int, 
                                  readloc::Int;
                                  ispos::Bool=true, 
                                  geneind::GeneInt=convert(GeneInt,searchsortedlast( lib.offset, indx )) )
   align = SGAlignment(zero(UInt32), zero(ReadLengthInt), SGAlignNode[], ispos, false)
   _ungapped_align!( p, lib, read, indx, readloc, align; ispos=ispos, geneind=geneind )
end

@inline function ungapped_align( p::AlignParam, 
                                 lib::SeqLibrary, 
                                 read::FASTQRecord;
                                 ispos::Bool=true, 
                                 anchor_left::Bool=true )

   seed,readloc,pos = seed_locate( p, lib.index, read, offset_left=anchor_left, both_strands=!p.is_stranded )
   isempty(seed) && return Nullable{Vector{SGAlignment}}()
   readlen = convert(ReadLengthInt, length(read.sequence))

   if !pos
      reverse_complement!( read )
      ispos = false
   end

   res      = Nullable{Vector{SGAlignment}}()
   maxscore = zero(Float32)

   for sidx in 1:length(seed)
      s = FMIndexes.sa_value( seed[sidx], lib.index ) + 1

      align = _ungapped_align( p, lib, read, s, readloc, ispos=ispos )
      scvar = identity(align, readlen)

      if isvalid(align)
         if isnull( res )
            res = Nullable(SGAlignment[ align ])
            maxscore = scvar
         else
            # new best score
            if scvar > maxscore
               # better than threshold for all of the previous scores
               if scvar - maxscore > p.score_range
                  length( res.value ) >= 1 && empty!( res.value )
               else # keep at least one previous score
                  splice_by_identity!( res.value, scvar, p.score_range, readlen )
               end
               maxscore = scvar
               push!( res.value, align )
            else # new score is lower than previously seen
               if maxscore - scvar <= p.score_range # but tolerable
                  push!( res.value, align )
               end
            end # end score vs maxscore
         end # end isnull
      end # end isvalid

   end
   res
end

@inline function spliced_fwd_extend( p::AlignParam,
                                     lib::SeqLibrary,
                                     read::FASTQRecord,
                                     geneind::CoordInt,
                                     edgeind::NodeInt,
                                     ridx::Int,
                                     right_kmer_ind::Int,
                                     align::SGAlignment )
   @inbounds begin
      edges_right = lib.edges.right
      isassigned( edges_right, right_kmer_ind ) || return align

      edges_left = lib.edges.left
      left_kmer_ind = kmer_index(lib.graphs[geneind].edgeleft[edgeind])
      isassigned( edges_left, left_kmer_ind ) || return align

      arrA = edges_left[ left_kmer_ind ]
      arrB = edges_right[ right_kmer_ind ]
      target_gene = align.path[1].gene
      rng = intersect_sorted_gene_range( arrA, arrB, target_gene )
      isempty(rng) && return align

      best = align
      best_score = score(best)
      init_score = best_score

      for idx in rng
         rn = arrB[idx]
         rn_offset = lib.graphs[rn.gene].nodeoffset[rn.node]

         res_align = ungapped_fwd_extend( p, lib, rn.gene, Int(rn_offset), read, ridx,
                                          align=copy_align(align), nodeidx=rn.node )
         res_score = score(res_align)
         if res_score > best_score
            best = res_align
            best_score = res_score
         end
      end
      if best_score >= init_score + p.kmer_size
         best.isvalid = true
      end
      return best
   end
end

@inline function spliced_rev_extend( p::AlignParam,
                                     lib::SeqLibrary,
                                     read::FASTQRecord,
                                     geneind::GeneInt,
                                     edgeind::NodeInt,
                                     ridx::Int,
                                     left_kmer_index::Int,
                                     align::SGAlignment )
   @inbounds begin
      edges_left = lib.edges.left
      isassigned( edges_left, left_kmer_index ) || return align

      edges_right = lib.edges.right
      right_kmer_ind = kmer_index(lib.graphs[geneind].edgeright[edgeind])
      isassigned( edges_right, right_kmer_ind ) || return align

      arrA = edges_right[ right_kmer_ind ]
      arrB = edges_left[ left_kmer_index ]
      target_gene = align.path[1].gene
      rng = intersect_sorted_gene_range( arrA, arrB, target_gene )
      isempty(rng) && return align

      best = align
      best_score = score(best)
      init_score = best_score

      for idx in rng
         rn = arrB[idx]
         rn_offset = lib.graphs[rn.gene].nodeoffset[rn.node-1] + lib.graphs[rn.gene].nodelen[rn.node-1] - 1

         res_align = ungapped_rev_extend( p, lib, rn.gene, Int(rn_offset), read, ridx,
                                          align=copy_align(align), nodeidx=rn.node-0x01 )
         res_score = score(res_align)
         if res_score > best_score
            best = res_align
            best_score = res_score
         end
      end
      if best_score >= init_score + p.kmer_size
         best.isvalid = true
      end
      return best
   end
end


# This is the main ungapped alignment extension function in the --> direction
# Returns: SGAlignment
@inbounds function ungapped_fwd_extend( p::AlignParam, 
                                        lib::SeqLibrary, 
                                        geneind::GeneInt, 
                                        sgidx::Int,
                                        read::FASTQRecord, 
                                        ridx::Int; 
                                        ispos::Bool=true,
                                        align::SGAlignment=SGAlignment(sgidx, ridx, SGAlignNode[],ispos,false),
                                        nodeidx::NodeInt=convert(NodeInt,searchsortedlast(lib.graphs[geneind].nodeoffset,CoordInt(sgidx))),
                                        init_mat::MatchesInt=zero(MatchesInt) )
   sg          = lib.graphs[geneind]
   rseq        = read.sequence
   rqual       = read.quality
   sgseq       = sg.seq
   sg_len      = length(sgseq)
   nodeoffset  = sg.nodeoffset
   nodelen     = sg.nodelen
   edgetype    = sg.edgetype
   readlen     = length(rseq)
   ksize       = p.kmer_size
   max_mis     = p.mismatches

   passed_edges::Union{Nothing, Vector{UInt8}} = nothing
   mistol      = mistolerance(align)

   push!( align.path, SGAlignNode( geneind, nodeidx, SGAlignScore(init_mat, zero(MismatchInt), zero(MismatchFloat)) ) ) # starting node
   cur_node    = nodeidx
   cur_mat     = init_mat
   cur_mis     = zero(MismatchInt)
   cur_mistol  = zero(MismatchFloat)
   need_flush  = true
   next_edge   = Int(nodeoffset[nodeidx]) + Int(nodelen[nodeidx])

   while mistol <= max_mis && ridx <= readlen && sgidx <= sg_len
      if sgidx == next_edge # hit next edge
         curedge = nodeidx + 0x01
         et = edgetype[curedge]
         if et == EDGETYPE_LR &&
            nodelen[curedge] >= ksize &&
            readlen - ridx + 1 >= ksize
            rkmer_ind = kmer_index_sub(rseq, ridx, ridx + ksize - 1)
            rkmer_ind == 0 && break
            if kmer_index(sg.edgeright[curedge]) == rkmer_ind
               align.path[end] = SGAlignNode( geneind, cur_node, SGAlignScore(cur_mat, cur_mis, cur_mistol) )
               ridx    += ksize - 1
               sgidx   += ksize - 1
               nodeidx += 0x01
               cur_node   = nodeidx
               cur_mat    = convert(MatchesInt, ksize - 1)
               cur_mis    = zero(MismatchInt)
               cur_mistol = zero(MismatchFloat)
               push!( align.path, SGAlignNode( geneind, nodeidx, SGAlignScore(cur_mat, cur_mis, cur_mistol) ) )
               align.isvalid = true
               next_edge = Int(nodeoffset[nodeidx]) + Int(nodelen[nodeidx])
            else
               align.path[end] = SGAlignNode( geneind, cur_node, SGAlignScore(cur_mat, cur_mis, cur_mistol) )
               need_flush = false
               align = spliced_fwd_extend( p, lib, read, geneind, curedge, ridx, rkmer_ind, align )
               break
            end
         elseif et == EDGETYPE_LR || et == EDGETYPE_LL
            align.path[end] = SGAlignNode( geneind, cur_node, SGAlignScore(cur_mat, cur_mis, cur_mistol) )
            if passed_edges === nothing
               passed_edges = Vector{UInt8}()
            end
            push!(passed_edges, UInt8(length(align.path)))
            nodeidx += 0x01
            cur_node   = nodeidx
            cur_mat    = zero(MatchesInt)
            cur_mis    = zero(MismatchInt)
            cur_mistol = zero(MismatchFloat)
            push!( align.path, SGAlignNode( geneind, nodeidx, zero(SGAlignScore) ) )
            next_edge = Int(nodeoffset[nodeidx]) + Int(nodelen[nodeidx])
         elseif et == EDGETYPE_LS
            if readlen - ridx + 1 >= ksize
               rkmer_ind = kmer_index_sub(rseq, ridx, ridx + ksize - 1)
               rkmer_ind == 0 && break
               align.path[end] = SGAlignNode( geneind, cur_node, SGAlignScore(cur_mat, cur_mis, cur_mistol) )
               need_flush = false
               align = spliced_fwd_extend( p, lib, read, geneind, curedge, ridx, rkmer_ind, align )
            end
            break
         elseif et == EDGETYPE_SR
            break
         else #'RR' || 'SL' || 'RS'
            align.path[end] = SGAlignNode( geneind, cur_node, SGAlignScore(cur_mat, cur_mis, cur_mistol) )
            nodeidx += 0x01
            if nodeidx <= length(nodelen)
               cur_node   = nodeidx
               cur_mat    = zero(MatchesInt)
               cur_mis    = zero(MismatchInt)
               cur_mistol = zero(MismatchFloat)
               push!( align.path, SGAlignNode( geneind, nodeidx, zero(SGAlignScore) ) )
               next_edge = Int(nodeoffset[nodeidx]) + Int(nodelen[nodeidx])
            end
         end
      end
      if eq_dna2_dna4(rseq, ridx, sgseq, sgidx)
         cur_mat += 0x01
      else
         prob = phred_to_prob( rqual[ridx] )
         @fastmath cur_mistol += prob
         @fastmath mistol += prob
         cur_mis += 0x01
      end
      ridx  += 1
      sgidx += 1
   end

   if need_flush
      align.path[end] = SGAlignNode( geneind, cur_node, SGAlignScore(cur_mat, cur_mis, cur_mistol) )
   end

   # if edgemat < K, spliced_extension for each in length(edges)
   if passed_edges !== nothing
      for c in length(passed_edges):-1:1  #most recent edge first
         cval = Int(passed_edges[c])
         if matches_range(align.path, cval + 1, length(align.path)) >= ksize + max_mis
            align.isvalid = true
            break
         else
            cur_edge = align.path[cval].node + one(NodeInt)
            resize!( align.path, cval ) # extension into current node failed

            cur_ridx = ridx - (sgidx - Int(nodeoffset[ cur_edge ]))
            (cur_ridx + ksize - 1) <= readlen || continue

            rkmer_ind = kmer_index_sub(rseq, cur_ridx, cur_ridx + ksize - 1)
            rkmer_ind == 0 && continue
            res_align = spliced_fwd_extend( p, lib, read, geneind, cur_edge,
                                            cur_ridx, rkmer_ind, align )
            if length(res_align.path) > length(align.path) # we added a valid node
               align = res_align > align ? res_align : align
            end
         end
      end
   end

   # clean up any bad trailing nodes
   while length(align.path) >= 1 && (align.path[end].meta.mismatches >= align.path[end].meta.matches ||
                                     Int(align.path[end].meta.matches) - Int(align.path[end].meta.mismatches) <= 1)
      pop!( align.path )
      length(align.path) == 0 && (align.isvalid = false)
   end

   if identity(align, readlen) >= p.score_min
      align.isvalid = true
   end
   align
end


# This is the main ungapped alignment extension function in the <-- direction
# Returns: SGAlignment
@inbounds function ungapped_rev_extend( p::AlignParam, 
                                        lib::SeqLibrary, 
                                        geneind::GeneInt, 
                                        sgidx::Int,
                                        read::FASTQRecord, 
                                        ridx::Int; 
                                        ispos::Bool=true,
                                        align::SGAlignment=SGAlignment(sgidx, ridx, SGAlignNode[], ispos, false),
                                        nodeidx::NodeInt=convert(NodeInt,searchsortedlast(lib.graphs[geneind].nodeoffset,CoordInt(sgidx))) )
   sg          = lib.graphs[geneind]
   rseq        = read.sequence
   rqual       = read.quality
   sgseq       = sg.seq
   nodeoffset  = sg.nodeoffset
   nodelen     = sg.nodelen
   edgetype    = sg.edgetype
   readlen     = length(rseq)
   ksize       = p.kmer_size
   max_mis     = p.mismatches

   passed_edges::Union{Nothing, Vector{UInt8}} = nothing
   mistol      = mistolerance(align)

   if length(align.path) <= 0 || align.path[1].gene != geneind || align.path[1].node != nodeidx
      pushfirst!( align.path, SGAlignNode( geneind, nodeidx, zero(SGAlignScore) ) ) # starting node if not already there
   end

   first_n     = align.path[1]
   cur_gene    = first_n.gene
   cur_node    = first_n.node
   cur_mat     = first_n.meta.matches
   cur_mis     = first_n.meta.mismatches
   cur_mistol  = first_n.meta.mistolerance
   need_flush  = true
   prev_edge   = Int(nodeoffset[nodeidx]) - 1

   while mistol <= max_mis && ridx > 0 && sgidx > 0
      if sgidx == prev_edge # hit right edge
         leftnode = nodeidx - 0x01
         et = edgetype[nodeidx]
         if et == EDGETYPE_LR &&
            nodelen[leftnode] >= ksize &&
            ridx >= ksize

            lkmer_ind = kmer_index_sub(rseq, ridx - ksize + 1, ridx)
            lkmer_ind == 0 && break
            if kmer_index(sg.edgeleft[nodeidx]) == lkmer_ind
               align.path[1] = SGAlignNode( cur_gene, cur_node, SGAlignScore(cur_mat, cur_mis, cur_mistol) )
               ridx    -= ksize - 1
               sgidx   -= ksize - 1
               nodeidx -= 0x01
               cur_gene   = geneind
               cur_node   = nodeidx
               cur_mat    = convert(MatchesInt, ksize - 1)
               cur_mis    = zero(MismatchInt)
               cur_mistol = zero(MismatchFloat)
               pushfirst!( align.path, SGAlignNode( geneind, nodeidx, SGAlignScore(cur_mat, cur_mis, cur_mistol) ) )
               align.isvalid = true
               prev_edge = Int(nodeoffset[nodeidx]) - 1
            else
               align.path[1] = SGAlignNode( cur_gene, cur_node, SGAlignScore(cur_mat, cur_mis, cur_mistol) )
               need_flush = false
               align = spliced_rev_extend( p, lib, read, geneind, nodeidx, ridx, lkmer_ind, align )
               break
            end
         elseif et == EDGETYPE_LR || et == EDGETYPE_RR
            align.path[1] = SGAlignNode( cur_gene, cur_node, SGAlignScore(cur_mat, cur_mis, cur_mistol) )
            if passed_edges === nothing
               passed_edges = Vector{UInt8}()
            end
            push!(passed_edges, UInt8(length(align.path) - 1))
            nodeidx -= 0x01
            cur_gene   = geneind
            cur_node   = nodeidx
            cur_mat    = zero(MatchesInt)
            cur_mis    = zero(MismatchInt)
            cur_mistol = zero(MismatchFloat)
            pushfirst!( align.path, SGAlignNode( geneind, nodeidx, zero(SGAlignScore) ) )
            prev_edge = Int(nodeoffset[nodeidx]) - 1
         elseif et == EDGETYPE_SR
            if ridx >= ksize
               lkmer_ind = kmer_index_sub(rseq, ridx - ksize + 1, ridx)
               lkmer_ind == 0 && break
               align.path[1] = SGAlignNode( cur_gene, cur_node, SGAlignScore(cur_mat, cur_mis, cur_mistol) )
               need_flush = false
               align = spliced_rev_extend( p, lib, read, geneind, nodeidx, ridx, lkmer_ind, align )
            end
            break
         elseif et == EDGETYPE_LS
            break
         else #'LL' || 'RS' || 'SL'
            align.path[1] = SGAlignNode( cur_gene, cur_node, SGAlignScore(cur_mat, cur_mis, cur_mistol) )
            nodeidx -= 0x01
            if nodeidx > 0
               cur_gene   = geneind
               cur_node   = nodeidx
               cur_mat    = zero(MatchesInt)
               cur_mis    = zero(MismatchInt)
               cur_mistol = zero(MismatchFloat)
               pushfirst!( align.path, SGAlignNode( geneind, nodeidx, zero(SGAlignScore) ) )
               prev_edge = Int(nodeoffset[nodeidx]) - 1
            end
         end
      end
      if eq_dna2_dna4(rseq, ridx, sgseq, sgidx)
         cur_mat += 0x01
      else
         prob = phred_to_prob( rqual[ridx] )
         @fastmath cur_mistol += prob
         @fastmath mistol += prob
         cur_mis += 0x01
      end
      ridx  -= 1
      sgidx -= 1
   end

   if need_flush
      align.path[1] = SGAlignNode( cur_gene, cur_node, SGAlignScore(cur_mat, cur_mis, cur_mistol) )
   end

   cur_ridx  = ridx + 1
   cur_sgidx = sgidx + 1
   # if passed_extend < K, spliced_extension for each in length(edges)
   if passed_edges !== nothing
      for c in length(passed_edges):-1:1  #most recent edge first
         cval = Int(passed_edges[c])
         if matches_range(align.path, 1, length(align.path) - (cval + 1)) >= ksize + max_mis
            align.isvalid = true
            break
         else
            cur_edge = align.path[end-cval].node
            deleteat!( align.path, 1:(length(align.path)-(cval+1)) ) # extension into current node failed
            cur_ridx  = ridx + (Int(nodeoffset[ cur_edge ]) - sgidx)
            cur_sgidx = Int(nodeoffset[ cur_edge ])
            (cur_ridx - ksize) > 0 || continue

            lkmer_ind = kmer_index_sub(rseq, cur_ridx - ksize, cur_ridx - 1)
            lkmer_ind == 0 && continue
            res_align = spliced_rev_extend( p, lib, read, geneind, cur_edge,
                                            cur_ridx - 1, lkmer_ind, align )
            if length(res_align.path) > length(align.path) # we added a valid node
               align = res_align > align ? res_align : align
            end
         end
      end
   end

   if identity(align, readlen) >= p.score_min
      align.isvalid = true
   end
   if sgidx < align.offset
      align.offset = convert(NodeInt, cur_sgidx)
   end
   if cur_ridx < align.offsetread
      align.offsetread = convert(ReadLengthInt, cur_ridx)
   end

   # clean up any bad left trailing nodes
   while length(align.path) >= 1 && (align.path[1].meta.mismatches >= align.path[1].meta.matches ||
                                     Int(align.path[1].meta.matches) - Int(align.path[1].meta.mismatches) <= 1)
      offset_change = align.path[1].meta.matches + align.path[1].meta.mismatches
      popfirst!( align.path )
      if length(align.path) >= 1
         align.offset      = nodeoffset[ align.path[1].node ]
         align.offsetread += offset_change
      else
         align.isvalid = false
      end
   end

   align
end
