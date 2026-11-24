
#=
#
#                     Basic IO Functions...
#
=#

plug_write( io::BufOut, str::S; plug::Char='\t' ) where {S <: AbstractString} = (write( io, str ); write( io, plug  ))
plug_write( io::BufOut, str::Char; plug::Char='\t' ) = (write( io, str ); write( io, plug  ))

tab_write( io::BufOut, str::S ) where {S <: AbstractString} = plug_write( io, str, plug='\t' )
tab_write( io::BufOut, str::Char ) = plug_write( io, str, plug='\t' )

end_write( io::BufOut, str::S ) where {S <: AbstractString} = plug_write( io, str, plug='\n' )
end_write( io::BufOut, str::Char ) = plug_write( io, str, plug='\n' )

function tab_write( io::BufOut, vs::Vector{SubString{String}} )
   for i in 1:length(vs)-1
      tab_write( io, vs[i] )
   end
   end_write( io, vs[length(vs)] )
end

const _FMT_BUF_POOL = [Vector{UInt8}(undef, 64) for _ in 1:64]

@inline function _fmt_buf()
   tid = Threads.threadid()
   if tid <= length(_FMT_BUF_POOL)
      return @inbounds _FMT_BUF_POOL[tid]
   else
      return Vector{UInt8}(undef, 64)
   end
end

@inline function write_int( io::BufOut, n::Integer )
   val = Int64(n)
   if val < 0
      write( io, UInt8('-') )
      u = reinterpret(UInt64, -val)
   else
      u = UInt64(val)
   end
   if u < 10
      write( io, UInt8('0') + UInt8(u) )
      return
   end
   buf = _fmt_buf()
   pos = 32
   @inbounds while u > 0
      q = div(u, UInt64(10))
      r = u - q * UInt64(10)
      buf[pos] = UInt8('0') + UInt8(r)
      pos -= 1
      u = q
   end
   GC.@preserve buf unsafe_write( io, pointer(buf, pos + 1), 32 - pos )
end

@inline function write_float( io::BufOut, x::Float64 )
   buf = _fmt_buf()
   pos = Base.Ryu.writeshortest( buf, 1, x )
   GC.@preserve buf unsafe_write( io, pointer(buf), pos - 1 )
end

function coord_write( io::BufOut, chr, first, last; tab=false )
   write( io, chr, UInt8(':') )
   write_int( io, first )
   write( io, UInt8('-') )
   write_int( io, last )
   tab && write( io, UInt8('\t') )
end

function coord_string( chr, first, last; tab::Bool=false )
   str = chr * ":" * string(first) * "-" * string(last)
   return tab ? str * "\t" : str
end

function coord_write( io::BufOut, chr, first, last, strand; tab=false )
   coord_write( io, chr, first, last )
   write( io, UInt8(':'), strand )
   tab && write( io, UInt8('\t') )
end

function complex_write( io::BufOut, complex::Int; tab=false )
   write( io, COMPLEX_CHAR )
   write_int( io, complex )
   tab && write( io, UInt8('\t') )
end

function conf_int_write( io::BufOut, conf_int::Tuple; tab=false, width=false, sig=4 )
   lo,hi = conf_int #unpack
   if width
      write_float( io, Float64(round(hi-lo, digits=sig)) )
      write( io, UInt8('\t') )
   end
   write_float( io, Float64(lo) )
   write( io, UInt8(',') )
   write_float( io, Float64(hi) )
   tab && write( io, UInt8('\t') )
end

@inline write_edge_node( io::BufOut, var::I ) where I <: Integer = write_int( io, var )
function write_edge_node( io::BufOut, var::NodeFloat )
   (root, node, pos) = decode_aberrant(var)
   write_int( io, root )
   root == node || write( io, UInt8('I') )
   if pos > 0
      write( io, UInt8('>') )
      write_int( io, pos - 1 )
   end
end

format_edge_string( var::I ) where I <: Integer = string(Int(var))

function format_edge_string( var::NodeFloat )
   (root, node, pos) = decode_aberrant(var)
   core_str = string(Int(root))
   core_str = root == node ? core_str : core_str * "I"
   if pos > 0
      return core_str * ">" * string(pos-1)
   else
      return core_str
   end
end

format_node_string( var::I ) where I <: Integer = string(Int(var))

function format_node_string( var::NodeFloat )
   (root, node, pos) = decode_aberrant(var)
   core_str = string(Int(root))
   core_str = root == node ? core_str : core_str * "I"
   if pos > 0
      return core_str * ">" * string(pos-1)
   else
      return core_str
   end
end

function edges_write( io::BufOut, edges::IntervalMap{N,C}, range ) where {N <: NodeNum, C <: ReadCounter}
   isempty(range) && return
   min_r = minimum(range)
   max_r = maximum(range)
   first_out = true
   for i in intersect( edges, (N(min_r), N(max_r)) )
      if i.first >= min_r && i.last <= max_r
         !first_out && write( io, UInt8(',') )
         write_edge_node( io, i.first )
         write( io, UInt8('-') )
         write_edge_node( io, i.last )
         write( io, UInt8(':') )
         write_float( io, Float64(round(get(i.value), digits=1)) )
         first_out = false
      end
   end
end

#=
#
#                     SAM Format IO...
#
=#

function seq_write( io::BufOut, read::FASTQRecord; tab=false )
   for c in read.sequence
      write( io, convert(Char, c) )
   end
   tab && write( io, '\t' )
end

function qual_write( io::BufOut, read::FASTQRecord; qualoffset=33, tab=false )
   for i in read.quality
      write( io, convert(Char, i+qualoffset) )
   end
   tab && write( io, '\t' )
end

function sam_flag( align::SGAlignment, lib::SeqLibrary, ind, paired, first, is_pair_rc, is_supplemental )
   flag = UInt16(0)
   if paired
      flag |= 0x01
      flag |= 0x02
      if first
         flag |= 0x40
         lib.info[ ind ].strand || (flag |= 0x10)
      else
         flag |= 0x80
         lib.info[ ind ].strand || (flag |= 0x10)
         #is_pair_rc   && (flag ⊻= 0x20)
      end
   else
      lib.info[ ind ].strand || (flag ⊻= 0x10)
   end
   is_supplemental && (flag |= 0x100)
   flag
end

@inline function soft_pad( align::SGAlignment )
   pad = align.offsetread - 1
   pad > 0 ? string(pad) * "S" : ""
end

function cigar_string( align::SGAlignment, sg::AnySpliceGraph, strand::Bool, readlen=0 )
   curpos    = align.offset                                 # left most position
   running   = 0                                            # cur number of running matches
   total     = align.offsetread - 1                         # total positions accounted for
   cigar     = String[ soft_pad(align) ]                    # build cigar string here

   for idx in 1:length(align.path)

      matches = align.path[idx].meta.matches + align.path[idx].meta.mismatches
      total  += matches
      curpos += matches

      if idx < length( align.path )
         curi = align.path[idx].node
         nexti = align.path[idx+1].node
         intron = strand ? Int(sg.nodecoord[nexti]) - Int(sg.nodecoord[curi] +  sg.nodelen[curi] - 1) - 1 :
                           Int(sg.nodecoord[curi])  - Int(sg.nodecoord[nexti] + sg.nodelen[nexti] - 1) - 1
         if intron >= 1
            push!( cigar, string(running + matches) * "M" )
            push!( cigar, string( intron ) * "N" )
            running = 0
            curpos  = sg.nodeoffset[nexti]
         else
            running += matches
         end
      else
         push!( cigar, string(running + matches) * "M" )
      end
   end
   if total < readlen
      soft = readlen - total
      push!( cigar, string( soft ) * "S" )
   end

   !strand && reverse!( cigar )
   join( cigar, "" ), curpos
end


# Write single SAM entry for one SGAlignment
function write_sam( io::BufOut, 
                  read::FASTQRecord, 
                  align::SGAlignment, 
                  lib::SeqLibrary;
                     mapq::Int=0, 
                     paired::Bool=false, 
                     fwd_mate::Bool=true, 
                     is_pair_rc::Bool=true, 
                     qualoffset::Int=33,
                     supplemental::Bool=false, 
                     tagstr::String="" )

   (align.isvalid && length(align.path) >= 1) || return
   geneind = align.path[1].gene
   strand  = lib.info[geneind].strand
   nodeind = strand ? align.path[1].node : align.path[end].node
   (align.path[end].node < nodeind || align.path[end].node < align.path[1].node) && return # TODO: allow circular SAM output
   sg = lib.graphs[geneind]
   cigar,endpos = cigar_string( align, sg, strand, length(read.sequence) )
   offset = strand ? (align.offset - sg.nodeoffset[nodeind]) : (sg.nodelen[nodeind] - (endpos - sg.nodeoffset[nodeind]))
   if !strand && !supplemental
      reverse_complement!(read)
   end
   tab_write( io, string(read) )
   tab_write( io, string( sam_flag(align, lib, geneind, paired, fwd_mate, is_pair_rc, supplemental) ) )
   tab_write( io, lib.info[geneind].name )
   tab_write( io, string( sg.nodecoord[nodeind] + offset ) )
   tab_write( io, string(mapq) )
   tab_write( io, cigar )
   tab_write( io, '*' )
   tab_write( io, '0' )
   tab_write( io, '0' )
   if !supplemental
      seq_write( io, read, tab=true )
      qual_write( io, read, qualoffset=qualoffset )
   else
      tab_write( io, '*' )
      write( io, '*' )
   end
   tagstr = "XS:A:" * (lib.info[geneind].strand ? '+' : '-')
   if tagstr != ""
      write( io, '\t' )
      write( io, tagstr )
   end
   write( io, '\n' )
end

# Efficient version for write_sam with multiple supplementary alignments
function indmax_score( vec::Vector{SGAlignment} )
   ind = 1
   max = score(vec[1])
   for i in 2:length(vec)
      cur = score(vec[i])
      if cur > max
         max = cur
         ind = i
      end
   end
   ind
end

function indmax_score( fwd::Vector{SGAlignment}, rev::Vector{SGAlignment} )
   ind = 1
   max = score(fwd[1], rev[1])
   for i in 2:length(fwd)
      cur = score(fwd[i], rev[i])
      if cur > max
         max = cur
         ind = i
      end
   end
   ind
end

# Write multiple SAM entries for a Vector of SGAlignment, highest score gets regular entry
# others get supplementary entries.
function write_sam( io::BufOut, 
                  read::FASTQRecord, 
                  alignvec::Vector{SGAlignment}, 
                  lib::SeqLibrary;
                     paired::Bool=false, 
                     fwd_mate::Bool=true, 
                     is_pair_rc::Bool=true, 
                     qualoffset::Int=33 )

   best = indmax_score( alignvec )
   mapq_val = 10 * length(alignvec)
   write_sam( io, read, alignvec[best], lib, mapq=mapq_val, paired=paired, fwd_mate=fwd_mate,
                        is_pair_rc=is_pair_rc, qualoffset=qualoffset, supplemental=false )
   for i in 1:length(alignvec)
      if i != best
         write_sam( io, read, alignvec[i], lib, mapq=mapq_val, paired=paired, fwd_mate=fwd_mate,
                        is_pair_rc=is_pair_rc, qualoffset=qualoffset, supplemental=true )

      end
   end
end

# Write multiple SAM entries for two vectors of SGAlignments for paired-end reads
function write_sam( io::BufOut, 
                  fwd::FASTQRecord, 
                  rev::FASTQRecord,
                  fwd_vec::Vector{SGAlignment}, 
                  rev_vec::Vector{SGAlignment}, 
                  lib::SeqLibrary;
                     paired::Bool=true, 
                     fwd_mate::Bool=true, 
                     is_pair_rc::Bool=true, 
                     qualoffset::Int=33 )

   best = indmax_score( fwd_vec, rev_vec )
   mapq_val = 10 * length(fwd_vec)
   write_sam( io, fwd, fwd_vec[best], lib,  mapq=mapq_val, paired=paired, fwd_mate=fwd_mate, is_pair_rc=is_pair_rc,
                                      qualoffset=qualoffset, supplemental=false )
   write_sam( io, rev, rev_vec[best], lib, mapq=mapq_val, paired=paired, fwd_mate=!fwd_mate, is_pair_rc=is_pair_rc,
                                      qualoffset=qualoffset, supplemental=false )
   for i in 1:length(fwd_vec)
      if i != best
         write_sam( io, fwd, fwd_vec[i], lib, mapq=mapq_val, paired=paired, fwd_mate=fwd_mate, is_pair_rc=is_pair_rc,
                                         qualoffset=qualoffset, supplemental=true )
         write_sam( io, rev, rev_vec[i], lib, mapq=mapq_val, paired=paired, fwd_mate=!fwd_mate, is_pair_rc=is_pair_rc,
                                         qualoffset=qualoffset, supplemental=true )
      end
   end
end


function write_sam_header( io::BufOut, lib::SeqLibrary, head::Bool=true )
   refs = Dict{String,Int}()
   for gind in 1:length(lib.graphs)
      name = lib.info[gind].name
      len  = Int(maximum( lib.graphs[gind].nodecoord ))
      if !haskey( refs, name ) || refs[name] < len
         refs[name] = len + 10000
      end
   end
   head && write( io, "@HD\tVN:1.0\tSO:unsorted\n" )
   for k in keys(refs)
      print( io, "@SQ\tSN:", k, "\tLN:", refs[k], '\n' )
   end
end

#=
#
#                     Event IO...
#
=#

function output_utr( io::BufOut, 
                     psi::Vector{Float64}, 
                     pgraph::Nullable{PsiGraph{S}},
                     total_reads::Float64, 
                     motif::String, 
                     sg::AnySpliceGraph, 
                     node::Int,
                     info::GeneMeta; 
                     empty=false,
                     empty_len::Int=length(psi) ) where {S <: MonotonicSet}

   st = motif == "TS" ? node : max(1, node - 1)
   en = st + (empty ? empty_len : length(psi)) - 1
   if en < st
      error("$psi, $total_reads, $motif, $node, $info, $empty")
   end
   i = 1
   @inbounds for n in st:en
      nc = sg.nodecoord[n]
      nl = sg.nodelen[n]
      write( io, info[1], UInt8('\t') )
      write_int( io, n )
      write( io, UInt8('\t'), info[2], UInt8(':') )
      write_int( io, nc )
      write( io, UInt8('-') )
      write_int( io, nc + nl - one(CoordInt) )
      write( io, UInt8('\t'), UInt8(info[3]), UInt8('\t'), motif, UInt8('\t') )
      if empty # had to add this flag since we iterate through the TS/TE nodes
         write( io, "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n" )
         i += 1
         continue
      end
      write_float( io, psi[i] )
      write( io, UInt8('\t') )
      if !isnull( pgraph )
         conf_int  = beta_posterior_ci( psi[i], total_reads, sig=3 )
         conf_int_write( io, conf_int, tab=true, width=true)
         write_float( io, round(total_reads, digits=3) )
         write( io, UInt8('\t') )
         complex_write( io, complexity( pgraph.value ), tab=true )
         write_float( io, round( shannon_index( pgraph.value ), digits=4 ) )
         write( io, UInt8('\t') )
         count_write( io, get(pgraph).nodes[i], get(pgraph).psi[i], tab=true )
         write( io, "NA\t" )
      else
         write( io, "NA\tNA\tNA\tNA\tNA\tNA\tNA\t\t" )
      end
      write( io, "NA\n" )
      i += 1
   end
   en
end

function output_psi( io::BufOut,
                     psi::Float64,
                     inc::Nullable{PsiGraph{S}}, 
                     exc::Nullable{PsiGraph{S}},
                     total_reads::Float64, 
                     conf_int::Tuple, 
                     motif::String,
                     sg::AnySpliceGraph,
                     edges::IntervalMap{N,C},
                     node::Int,
                     info::GeneMeta,
                     bias ) where {N <: NodeNum, C <: ReadCounter, S <: MonotonicSet}
   @inbounds begin
      nlen = sg.nodelen[node]
      nlen == 0 && return
      left  = sg.nodecoord[node]
      right = left + nlen - one(NodeInt)
   end
   write( io, info[1], UInt8('\t') )
   write_int( io, node )
   write( io, UInt8('\t'), info[2], UInt8(':') )
   write_int( io, left )
   write( io, UInt8('-') )
   write_int( io, right )
   write( io, UInt8('\t'), UInt8(info[3]), UInt8('\t'), motif, UInt8('\t') )
   write_float( io, psi )
   write( io, UInt8('\t') )
   _output_psi_tail( io, psi, inc, exc, total_reads, conf_int, edges )
end

function output_psi( io::BufOut, 
                     psi::Float64, 
                     inc::Nullable{PsiGraph{S}}, 
                     exc::Nullable{PsiGraph{S}},
                     total_reads::Float64, 
                     conf_int::Tuple, 
                     motif::String,
                     coord::String,
                     edges::IntervalMap{N,C},
                     node, 
                     info::GeneMeta, 
                     bias ) where {N <: NodeNum, C <: ReadCounter, S <: MonotonicSet}

   write( io, info[1], UInt8('\t') )
   write_edge_node( io, node )
   write( io, UInt8('\t'), coord, UInt8('\t'), UInt8(info[3]), UInt8('\t'), motif, UInt8('\t') )
   write_float( io, psi )
   write( io, UInt8('\t') )
   _output_psi_tail( io, psi, inc, exc, total_reads, conf_int, edges )
end

function _output_psi_tail( io::BufOut,
                           psi::Float64,
                           inc::Nullable{PsiGraph{S}},
                           exc::Nullable{PsiGraph{S}},
                           total_reads::Float64,
                           conf_int::Tuple,
                           edges::IntervalMap{N,C} ) where {N <: NodeNum, C <: ReadCounter, S <: MonotonicSet}
   conf_int_write( io, conf_int, tab=true, width=true )
   write_float( io, round(total_reads, digits=3) )
   write( io, UInt8('\t') )

   if !isnull( inc ) && !isnull( exc )
      complex_write( io, complexity( inc.value, exc.value ), tab=true )
      write_float( io, round( shannon_index( inc.value, exc.value ), digits=4 ) )
      write( io, UInt8('\t') )
   else
      write( io, "K0\t0.0\t" )
   end

   min_n,max_n = 10000000.0,0.0
   if !isnull( inc )
      count_write( io, get(inc), psi, tab=true )
      min_n,max_n = get(inc).min,get(inc).max
   else
      write( io, "NA\t" )
   end

   if !isnull( exc )
      count_write( io, get(exc), 1.0-psi, tab=true )
      min_n = get(exc).min < min_n ? get(exc).min : min_n
      max_n = get(exc).max > max_n ? get(exc).max : max_n
   else
      write( io, "NA\t" )
   end

   edges_write( io, edges, min_n:max_n )
   write( io, UInt8('\n') )
end

function output_circular( io::BufOut, sg::AnySpliceGraph, sgquant::SpliceGraphQuant, info::GeneMeta )
   for (st,en) in keys(sgquant.circ)
      back_len = 0.0
      back_cnt = 0.0
      fore_cnt = get(sgquant.circ[(st,en)])
      for edg in intersect( sgquant.edge, (st, st) )
         if edg.first == st
            back_len += 1
            back_cnt += get(edg.value)
         end
      end
      total_reads = (fore_cnt + back_cnt)
      psi = fore_cnt / total_reads
      if !isnan(psi) && total_reads > 0
         write( io, info[1], UInt8('\t') )
         write_int( io, st )
         write( io, UInt8('-') )
         write_int( io, en )
         write( io, UInt8('\t') )
         coord_write( io, info[2], sg.nodecoord[st]+sg.nodelen[st]-1, sg.nodecoord[en], tab=true )
         write( io, UInt8(info[3]), "\tBS\t" )
         write_float( io, psi )
         write( io, UInt8('\t') )
         conf_int  = beta_posterior_ci( psi, total_reads, sig=3 )
         conf_int_write( io, conf_int, tab=true, width=true)
         write_float( io, round(total_reads, digits=3) )
         write( io, "\tNA\tNA\tNA\tNA\tNA\n" )
      end
   end
end



function output_empty( io::BufOut, motif::String, sg::AnySpliceGraph, node::Int, info::GeneMeta )
   @inbounds begin
      nlen = sg.nodelen[node]
      nlen == 0 && return
      ncoord = sg.nodecoord[node]
      write( io, info[1], UInt8('\t') )
      write_int( io, node )
      write( io, UInt8('\t'), info[2], UInt8(':') )
      write_int( io, ncoord )
      write( io, UInt8('-') )
      write_int( io, ncoord + nlen - one(CoordInt) )
      write( io, UInt8('\t'), UInt8(info[3]), UInt8('\t'), motif, "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n" )
   end
end

function write_monotonic_set( io::BufOut, is::S ) where S <: MonotonicSet
   first_el = true
   for n in is
      !first_el && write( io, UInt8('-') )
      write_edge_node( io, n )
      first_el = false
   end
end

function Base.string( is::S ) where S <: MonotonicSet
   buf = IOBuffer()
   first_el = true
   for n in is
      !first_el && write( buf, '-' )
      write( buf, format_edge_string(n) )
      first_el = false
   end
   String(take!(buf))
end

function count_write( io::BufOut, nodestr::S, psi::Float64; tab=false ) where S <: MonotonicSet
   write_monotonic_set( io, nodestr )
   write( io, UInt8(':') )
   write_float( io, round(psi, digits=6) )
   tab && write( io, UInt8('\t') )
end

function count_write( io::BufOut, nodestr, psi::Float64; tab=false )
   print( io, nodestr )
   write( io, UInt8(':') )
   write_float( io, round(psi, digits=6) )
   tab && write( io, UInt8('\t') )
end

function count_write( io::BufOut, nodestr, psi::String; tab=false )
   print( io, nodestr )
   write( io, UInt8(':'), psi )
   tab && write( io, UInt8('\t') )
end

function count_write( io::BufOut, pgraph::PsiGraph, psi::Float64; tab=false )
   start = 1
   n_nodes = length(pgraph.nodes)
   if psi == 1.0 && n_nodes > 1
      @inbounds for i in 1:n_nodes-1
         write_monotonic_set( io, pgraph.nodes[i] )
         write( io, UInt8(',') )
      end
      start = n_nodes
   end
   @inbounds for i in start:n_nodes
      curpsi = i > length(pgraph.psi) ? psi : pgraph.psi[i]
      count_write( io, pgraph.nodes[i], curpsi )
      (i < n_nodes) && write( io, UInt8(',') )
   end
   tab && write( io, UInt8('\t') )
end

function output_psi_header( io::BufOut )
   tab_write( io, "Gene\tNode\tCoord\tStrand" )
   tab_write( io, "Type\tPsi\tCI_Width\tCI_Lo,Hi" )
   write( io, "Total_Reads\tComplexity\tEntropy\tInc_Paths\tExc_Paths\tEdges\n" )
end

function output_tpm_header( io::BufOut, idstr::String="Gene" )
   write( io, idstr )
   write( io, "\tTpM\tRead_Counts\n" )
end

function output_diff_header( io::BufOut )
   tab_write( io, "Gene\tNode\tCoord\tStrand" )
   write( io, "Type\tComplexity\tEntropy\tSampleCnt_A\tSampleCnt_B\tMeanReads_A\tMeanReads_B\tPsi_A\tPsi_B\tDeltaPsi\tProbability\n" )
end

function output_diff( io::BufOut, 
                     event, 
                     complex::Int, 
                     entropy::Float64,
		               samp_a::Int,
                     samp_b::Int,
                     mean_a::Float64,
                     mean_b::Float64,
                     psi_a::Float64, 
                     psi_b::Float64, 
                     deltapsi::Float64, 
                     prob::Float64, 
                     sig=5 )

   for i in 1:length(event)
      tab_write( io, event[i] )
   end
   write( io, COMPLEX_CHAR )
   tab_write( io, string(complex))
   tab_write( io, string(entropy) )
   tab_write( io, string(samp_a))
   tab_write( io, string(samp_b))
   tab_write( io, string(round(mean_a, digits=1)) )
   tab_write( io, string(round(mean_b, digits=1)) )
   tab_write( io, string(round(psi_a, digits=sig)) )
   tab_write( io, string(round(psi_b, digits=sig)) )
   tab_write( io, string(round(deltapsi, digits=sig)) )
   end_write( io, string(round(prob, digits=10)) )
end

#=
#
#                     Mapping Stats IO...
#
=#

function count_annotated_edges( sg::AnySpliceGraph, sgquant::SpliceGraphQuant )
   annoreads  = 0.0
   annocnt    = 0
   totalreads = 0.0
   totalcnt   = 0
   apaths = sg.annopath
   for edg in sgquant.edge
      v = get(edg.value)
      totalreads += v
      totalcnt   += 1
      if edg in apaths
         annoreads += v
         annocnt   += 1
      end
   end
   annoreads,annocnt,totalreads,totalcnt
end

function count_read_types( lib::SeqLibrary, graphq::GraphLibQuant )
   bodyreads = 0.0
   juncreads = 0.0
   annoreads = 0.0
   junccnt   = 0
   annocnt   = 0
   @inbounds for i in 1:length(lib.graphs)
      sgq = graphq.quant[i]
      for node in sgq.node
         bodyreads += get(node)
      end
      if !isempty(sgq.edge)
         areads,acnt,treads,tcnt = count_annotated_edges( lib.graphs[i], sgq )
         annoreads += areads
         juncreads += treads
         annocnt   += acnt
         junccnt   += tcnt
      end
   end
   Int(round(bodyreads)),Int(round(juncreads)),Int(round(annoreads)),junccnt,annocnt
end

function output_stats( filename::String, 
                       lib::SeqLibrary, 
                       graphq::GraphLibQuant, 
                       param::AlignParam,
                       index::String, 
                       total::Int, 
                       mapped::Int, 
                       multi::Int, 
                       readlen::Int, 
                       ver::String )

   io = open( filename, "w" )
   stream = ZlibDeflateOutputStream( io )

   output_stats( stream, lib, graphq, param, index, total, mapped, multi, readlen, ver )

   close(stream)
   close(io)
end

function output_stats( io::BufOut, 
                       lib::SeqLibrary, 
                       graphq::GraphLibQuant, 
                       param::AlignParam,
                       index::String, 
                       total::Int, 
                       mapped::Int, 
                       multi::Int, 
                       readlen::Int, 
                       ver::String )

   write( io, "# -------- Whippet $ver ---------\n" )
   write( io, "Mapped_Index\t$index\n" )
   write( io, "Read_Length\t$readlen\n" )
   write( io, "Total_Reads\t$total\n" )
   write( io, "# -------- Mapping Stats ---------\n" )
   write( io, "Mapped_Reads\t$mapped\n" )
   write( io, "Mapped_Percent\t$(round((mapped/total)*100,digits=2))%\n" )
   write( io, "# -------- Multi-loci Mapping ---------\n" )
   write( io, "Multimap_Reads\t$multi\n" )
   write( io, "Multimap_Percent\t$(round((multi/mapped)*100,digits=2))%\n" )
   write( io, "# -------- Unspliced Reads ---------\n" )
   body,junc,anno,jcnt,acnt = count_read_types( lib, graphq )
   write( io, "Exon_Body_Reads\t$body\n" )
   write( io, "Exon_Body_Percent\t$(round((body/(body+junc))*100,digits=2))%\n" )
   write( io, "# -------- Spliced Reads ---------\n" )
   write( io, "Junction_Reads\t$junc\n" )
   write( io, "Junction_Percent\t$(round((junc/(body+junc))*100,digits=2))%\n" )
   write( io, "Junction_Number\t$jcnt\n" )
   write( io, "Annotated_Junc_Reads\t$anno\n" )
   write( io, "Annotated_Junc_Percent\t$(round((anno/junc)*100,digits=2))%\n" )
   write( io, "Annotated_Junc_Number\t$acnt\n" )
   write( io, "Novel_Junc_Number\t$(jcnt-acnt)\n" )
   write( io, "Novel_Junc_Percent\t$(round(((jcnt-acnt)/jcnt)*100,digits=2))%\n" )
end


function output_junctions( filename::String, lib::SeqLibrary, graphq::GraphLibQuant )
   io = open( filename, "w" )
   stream = ZlibDeflateOutputStream( io )

   output_junctions( stream, lib, graphq )

   close(stream)
   close(io)
end

function output_junctions( io::BufOut, lib::SeqLibrary, graphq::GraphLibQuant )
   function write_junction( sg::AnySpliceGraph, edg::IntervalValue, i::Int, str::String )
      info_i = @inbounds lib.info[i]
      if info_i.strand
         donor    = @inbounds sg.nodecoord[edg.first]+sg.nodelen[edg.first]-1
         acceptor = @inbounds sg.nodecoord[edg.last]
      else
         donor    = @inbounds sg.nodecoord[edg.last]+sg.nodelen[edg.last]-1
         acceptor = @inbounds sg.nodecoord[edg.first]
      end
      write( io, info_i.name, UInt8('\t') )
      write_int( io, donor )
      write( io, UInt8('\t') )
      write_int( io, acceptor )
      write( io, UInt8('\t'), info_i.gene, UInt8(':') )
      write_int( io, edg.first )
      write( io, UInt8('-') )
      write_int( io, edg.last )
      write( io, UInt8(':'), str, UInt8('\t') )
      write_float( io, get(edg.value) )
      write( io, UInt8('\t'), info_i.strand ? UInt8('+') : UInt8('-'), UInt8('\n') )
   end

   function write_junctions( sg::AnySpliceGraph, sgquant::SpliceGraphQuant, i::Int )
      isempty(sgquant.edge) && return
      et = sg.edgetype
      ap = sg.annopath
      @inbounds for edg in sgquant.edge
         if et[edg.first+1] in (EDGETYPE_LS, EDGETYPE_LL, EDGETYPE_LR) &&
            et[edg.last]    in (EDGETYPE_SR, EDGETYPE_RR, EDGETYPE_LR)

            if inall( edg, ap )
               write_junction( sg, edg, i, "CON_ANNO" )
            elseif edg in ap
               write_junction( sg, edg, i, "ALT_ANNO" )
            else # unique edge
               write_junction( sg, edg, i, "ALT_UNIQ" )
            end
         end
      end
   end

   write( io, "CHROM\tDONOR_COORD\tACCEPTOR_COORD\tGENEID:NODES:TYPE\tREADCOUNT\tSTRAND\n")
   @inbounds for i in 1:length(lib.graphs)
      write_junctions( lib.graphs[i], graphq.quant[i], i )
   end
end


function output_exons( nodecoord::Vector{CoordInt},
                       nodelen::Vector{CoordInt},
                       tx::RefTx, 
                       strand::Bool )
   path = BitSet()
   # this may be `poor form`, but 256 is too big for default!
   # resize!(path.bits, 64) # Deprecated:  = zeros(UInt32,64>>>5)
   for i in 1:length(tx.acc)
      ind = collect(searchsorted( nodecoord, tx.acc[i], rev=!strand ))[end]
      push!( path, ind )
      cur = ind + (strand ? 1 : -1)
      while 1 <= cur <= length(nodecoord) &&
            nodecoord[cur] <= tx.don[i]
         push!( path, cur )
         cur = cur + (strand ? 1 : -1 )
      end
   end
   path
end

function Base.in( nodes::UnitRange, is::S ) where S <: MonotonicSet
   for i in nodes
      if !(i in is)
         return false
      end
   end
   return true
end

function Base.in( nodes::UnitRange, v::Vector{S} ) where S <: MonotonicSet
   for is in v
      if nodes in is
         return true
      end
   end
   return false
end

function output_exons( io::BufOut, sg::SpliceGraph, info::GeneInfo )
   for i in 1:length(sg.nodecoord)
      if isexonstart(sg.edgetype[i])
         for j in i:length(sg.nodecoord)
            if isexonstop(sg.edgetype[j+1])
               # print exon
               anno = i:j in sg.annopath ? "Y" : "N"
               tab_write( io, info.gene )
               if info.strand
                  coord_write( io, info.name, sg.nodecoord[i], sg.nodecoord[j] + sg.nodelen[j] - 1, info.strand ? "+" : "-", tab=true )
               else
                  coord_write( io, info.name, sg.nodecoord[j], sg.nodecoord[i] + sg.nodelen[i] - 1, info.strand ? "+" : "-", tab=true )
               end
               tab_write( io, join(i:j, ",") )
               end_write( io, anno )
            end
            if ishardedge(sg.edgetype[j+1])
               break
            end
         end
      end
   end
end

function output_exons( io::BufOut, lib::GraphLib )
   for i in 1:length(lib.graphs)
      output_exons( io, lib.graphs[i], lib.info[i] )
   end
end

function output_exons( filename::String, lib::GraphLib )
   io = open( filename, "w" )
   stream = ZlibDeflateOutputStream( io )

   write( stream, "Gene\tPotential_Exon\tWhippet_Nodes\tIs_Annotated\n" )
   output_exons( stream, lib )

   close(stream)
   close(io)
end

function output_paths( io::BufOut, sg::SpliceGraph, info::GeneInfo )
   for i in 1:length(sg.annopath)
      tab_write( io, info.gene )
      tab_write( io, sg.annoname[i] )
      tab_write( io, join(sg.annopath[i], '-') )

      next = iterate(sg.annopath[i])
      prev = 0
      while next !== nothing
         (val, state) = next
         next = iterate(sg.annopath[i], state)
         if prev != 0
            write( io, string(prev) )
            write( io, '-' )
            write( io, string(val) )
            next !== nothing && write(io, ',')
         end
         prev = val
      end
      write( io, '\n' )
   end
end

function output_paths( io::BufOut, lib::GraphLib )
   for i in 1:length(lib.graphs)
      output_paths( io, lib.graphs[i], lib.info[i] )
   end
end

function output_paths( filename::String, lib::GraphLib )
   io = open( filename, "w" )
   stream = ZlibDeflateOutputStream( io )

   write( stream, "Gene\tIsoform\tPath\tEdges\n" )
   output_paths( stream, lib )

   close(stream)
   close(io)
end
