
@inline score( one::A, two::A  ) where {A <: UngappedAlignment} = @fastmath score( one ) + score( two )
@inline identity( one::A, two::A, onelen::Int, twolen::Int ) where {A <: UngappedAlignment} = @fastmath score( one, two ) / ( onelen + twolen )

@inline function sorted_offsets!( res::Vector{Int}, sa::UnitRange, index::FMIndexes.FMIndex )
   n = length(sa)
   resize!(res, n)
   ridx = 1
   @inbounds for sidx in sa
      res[ridx] = FMIndexes.sa_value( sidx, index ) + 1
      ridx += 1
   end
   n > 1 && sort!(res)
   res
end

@inline function sorted_offsets( sa::UnitRange, index::FMIndexes.FMIndex )
   sorted_offsets!( Vector{Int}(undef, length(sa)), sa, index )
end

# Paired version.
function splice_by_score!( one::Vector{A}, two::Vector{A}, threshold, buffer ) where A <: UngappedAlignment
   i = 1
   while i <= length( one )
      if threshold - score( one[i], two[i] ) > buffer
         splice!( one, i )
         splice!( two, i )
         i -= 1
      end
      i += 1
   end
end

function splice_by_identity!( one::Vector{A}, two::Vector{A}, threshold, buffer, onelen, twolen ) where A <: UngappedAlignment
   i = 1
   while i <= length( one )
      if threshold - identity( one[i], two[i], onelen, twolen ) > buffer
         splice!( one, i )
         splice!( two, i )
         i -= 1
      end
      i += 1
   end
end

# Paired ungapped_align
function ungapped_align( p::AlignParam, lib::SeqLibrary, fwd::FASTQRecord, rev::FASTQRecord;
                         ispos::Bool=true, anchor_left::Bool=true,
                         fwd_buf::Vector{Int}=Int[], rev_buf::Vector{Int}=Int[] )

   empty_res = Nullable{Vector{SGAlignment}}()
   fwd_seed,fwd_readloc,strand = seed_locate( p, lib.index, fwd, offset_left=anchor_left, both_strands=!p.is_stranded )
   isempty(fwd_seed) && return empty_res, empty_res

   if (p.is_pair_rc && strand) || (!p.is_pair_rc && !strand)
      reverse_complement!(rev)
      rev_seed,rev_readloc = seed_locate( p, lib.index, rev, offset_left=!anchor_left, both_strands=false)
   else
      rev_seed,rev_readloc = seed_locate( p, lib.index, rev, offset_left=anchor_left, both_strands=false)
   end
   isempty(rev_seed) && return empty_res, empty_res

   fwd_len = length(fwd.sequence)
   rev_len = length(rev.sequence)

   if !strand
      reverse_complement!(fwd)
      ispos = false
   end

   # Fast path: uniquely seeded pair (most common case, zero buffer allocation / sorting)
   if length(fwd_seed) == 1 && length(rev_seed) == 1
      fwd_off = FMIndexes.sa_value( first(fwd_seed), lib.index ) + 1
      rev_off = FMIndexes.sa_value( first(rev_seed), lib.index ) + 1
      if abs( fwd_off - rev_off ) < p.seed_pair_range
         geneind = convert( CoordInt, searchsortedlast( lib.offset, fwd_off ) )
         next_offset = geneind < length(lib.offset) ? (@inbounds lib.offset[geneind+1]) :
                                                      (@inbounds lib.offset[geneind] + length(lib.graphs[geneind].seq))
         if (@inbounds lib.offset[geneind]) <= rev_off < next_offset
            fwd_aln = _ungapped_align( p, lib, fwd, fwd_off, fwd_readloc; ispos=ispos, geneind=geneind )
            rev_aln = _ungapped_align( p, lib, rev, rev_off, rev_readloc; ispos=!ispos, geneind=geneind )
            if isvalid(fwd_aln) && isvalid(rev_aln)
               return Nullable(SGAlignment[ fwd_aln ]), Nullable(SGAlignment[ rev_aln ])
            end
         end
      end
      return empty_res, empty_res
   end

   fwd_res  = empty_res
   rev_res  = empty_res
   maxscore = 0.0f0

   fwd_sorted = sorted_offsets!( fwd_buf, fwd_seed, lib.index )
   rev_sorted = sorted_offsets!( rev_buf, rev_seed, lib.index )

   fidx,ridx = 1,1
   n_fwd = length(fwd_sorted)
   n_rev = length(rev_sorted)

   @inbounds while fidx <= n_fwd && ridx <= n_rev

      foff = fwd_sorted[fidx]
      roff = rev_sorted[ridx]
      dist = abs( foff - roff )
      if dist < p.seed_pair_range # inside allowable pair distance

         geneind = convert( CoordInt, searchsortedlast( lib.offset, foff ) )
         next_offset = geneind < length(lib.offset) ? lib.offset[geneind+1] :
                                                           lib.offset[geneind] + length(lib.graphs[geneind].seq)
         if lib.offset[geneind] <= roff < next_offset

            fwd_aln = _ungapped_align( p, lib, fwd, foff, fwd_readloc; ispos=ispos, geneind=geneind )
            rev_aln = _ungapped_align( p, lib, rev, roff, rev_readloc; ispos=!ispos, geneind=geneind )

            @fastmath scvar = identity( fwd_aln, rev_aln, fwd_len, rev_len )

            if isvalid(fwd_aln) && isvalid(rev_aln)
               if isnull( fwd_res ) || isnull( rev_res )
                  fwd_res = Nullable(SGAlignment[ fwd_aln ])
                  rev_res = Nullable(SGAlignment[ rev_aln ])
                  maxscore = scvar
               else
               # new best score
                  if scvar > maxscore
                     # better than threshold for all of the previous scores
                     if scvar - maxscore > p.score_range
                        if length( fwd_res.value ) >= 1
                           empty!( fwd_res.value )
                           empty!( rev_res.value )
                        end
                     else # keep at least one previous score
                        splice_by_identity!( fwd_res.value, rev_res.value, scvar, p.score_range, fwd_len, rev_len )
                     end
                     maxscore = scvar
                     push!( fwd_res.value, fwd_aln )
                     push!( rev_res.value, rev_aln )
                  else # new score is lower than previously seen
                     if maxscore - scvar <= p.score_range # but tolerable
                        push!( fwd_res.value, fwd_aln )
                        push!( rev_res.value, rev_aln )
                     end
                  end # end score vs maxscore
               end # end isnull
            end # end isvalid
            ridx += 1
            fidx += 1
            continue
         end # is reverse in same gene as fwd
      end # end dist vs. seed_pair_range
      if foff > roff
         ridx += 1
      else
         fidx += 1
      end
   end #end while

   fwd_res,rev_res
end

# Zero-allocation paired alignment writing directly into pre-allocated chunk slots.
# Returns:
#   Int32(0) -> unmapped
#   Int32(1) -> unique mapping stored in fwd_alns[slot] and rev_alns[slot]
#   Int32(k) -> k > 1 multi-mappings stored in fwd_multi[1:k] and rev_multi[1:k]
function ungapped_align_into!( p::AlignParam,
                               lib::SeqLibrary,
                               fwd::FASTQRecord,
                               rev::FASTQRecord,
                               fwd_alns::Vector{SGAlignment},
                               rev_alns::Vector{SGAlignment},
                               fwd_multi::Vector{SGAlignment},
                               rev_multi::Vector{SGAlignment},
                               slot::Int;
                               ispos::Bool=true,
                               anchor_left::Bool=true,
                               fwd_buf::Vector{Int}=Int[],
                               rev_buf::Vector{Int}=Int[] )::Int32

   fwd_seed,fwd_readloc,strand = seed_locate( p, lib.index, fwd, offset_left=anchor_left, both_strands=!p.is_stranded )
   isempty(fwd_seed) && return Int32(0)

   if (p.is_pair_rc && strand) || (!p.is_pair_rc && !strand)
      reverse_complement!(rev)
      rev_seed,rev_readloc = seed_locate( p, lib.index, rev, offset_left=!anchor_left, both_strands=false )
   else
      rev_seed,rev_readloc = seed_locate( p, lib.index, rev, offset_left=anchor_left, both_strands=false )
   end
   isempty(rev_seed) && return Int32(0)

   fwd_len = length(fwd.sequence)
   rev_len = length(rev.sequence)

   if !strand
      reverse_complement!(fwd)
      ispos = false
   end

   # Fast path: uniquely seeded pair (zero heap allocations in the non-branching case)
   if length(fwd_seed) == 1 && length(rev_seed) == 1
      fwd_off = FMIndexes.sa_value( first(fwd_seed), lib.index ) + 1
      rev_off = FMIndexes.sa_value( first(rev_seed), lib.index ) + 1
      if abs( fwd_off - rev_off ) < p.seed_pair_range
         geneind = convert( CoordInt, searchsortedlast( lib.offset, fwd_off ) )
         next_offset = geneind < length(lib.offset) ? (@inbounds lib.offset[geneind+1]) :
                                                      (@inbounds lib.offset[geneind] + length(lib.graphs[geneind].seq))
         if (@inbounds lib.offset[geneind]) <= rev_off < next_offset
            fa = _ungapped_align!( p, lib, fwd, fwd_off, fwd_readloc, (@inbounds fwd_alns[slot]); ispos=ispos, geneind=geneind )
            ra = _ungapped_align!( p, lib, rev, rev_off, rev_readloc, (@inbounds rev_alns[slot]); ispos=!ispos, geneind=geneind )
            @inbounds fwd_alns[slot] = fa
            @inbounds rev_alns[slot] = ra
            if isvalid(fa) && isvalid(ra)
               return Int32(1)
            end
         end
      end
      return Int32(0)
   end

   empty!(fwd_multi)
   empty!(rev_multi)
   maxscore = 0.0f0

   fwd_sorted = sorted_offsets!( fwd_buf, fwd_seed, lib.index )
   rev_sorted = sorted_offsets!( rev_buf, rev_seed, lib.index )

   fidx,ridx = 1,1
   n_fwd = length(fwd_sorted)
   n_rev = length(rev_sorted)

   @inbounds while fidx <= n_fwd && ridx <= n_rev
      foff = fwd_sorted[fidx]
      roff = rev_sorted[ridx]
      dist = abs( foff - roff )
      if dist < p.seed_pair_range
         geneind = convert( CoordInt, searchsortedlast( lib.offset, foff ) )
         next_offset = geneind < length(lib.offset) ? lib.offset[geneind+1] :
                                                      lib.offset[geneind] + length(lib.graphs[geneind].seq)
         if lib.offset[geneind] <= roff < next_offset
            fwd_aln = _ungapped_align( p, lib, fwd, foff, fwd_readloc; ispos=ispos, geneind=geneind )
            rev_aln = _ungapped_align( p, lib, rev, roff, rev_readloc; ispos=!ispos, geneind=geneind )
            @fastmath scvar = identity( fwd_aln, rev_aln, fwd_len, rev_len )

            if isvalid(fwd_aln) && isvalid(rev_aln)
               if isempty(fwd_multi)
                  push!(fwd_multi, fwd_aln)
                  push!(rev_multi, rev_aln)
                  maxscore = scvar
               else
                  if scvar > maxscore
                     if scvar - maxscore > p.score_range
                        empty!(fwd_multi)
                        empty!(rev_multi)
                     else
                        splice_by_identity!( fwd_multi, rev_multi, scvar, p.score_range, fwd_len, rev_len )
                     end
                     maxscore = scvar
                     push!(fwd_multi, fwd_aln)
                     push!(rev_multi, rev_aln)
                  else
                     if maxscore - scvar <= p.score_range
                        push!(fwd_multi, fwd_aln)
                        push!(rev_multi, rev_aln)
                     end
                  end
               end
            end
            ridx += 1
            fidx += 1
            continue
         end
      end
      if foff > roff
         ridx += 1
      else
         fidx += 1
      end
   end

   nh = Int32(length(fwd_multi))
   if nh == Int32(1)
      @inbounds fwd_alns[slot] = fwd_multi[1]
      @inbounds rev_alns[slot] = rev_multi[1]
      empty!(fwd_multi)
      empty!(rev_multi)
   end
   return nh
end
