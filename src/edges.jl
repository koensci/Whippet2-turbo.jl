
const SGNodeSet = Vector{SGNode}

struct Edges{K}
   left::Vector{SGNodeSet}
   right::Vector{SGNodeSet}
end

Base.:<( a::SGNode, b::SGNode ) = <( a.gene, b.gene )
Base.:>( a::SGNode, b::SGNode ) = >( a.gene, b.gene )
Base.:<=( a::SGNode, b::SGNode ) = <=( a.gene, b.gene )
Base.:>=( a::SGNode, b::SGNode ) = >=( a.gene, b.gene )

sortlt( a::SGNode, b::SGNode ) = a.gene == b.gene ? <( a.node, b.node ) : <( a.gene, b.gene )

Base.convert( ::Type{Edges{K}}, graphs::Vector{SpliceGraph} ) where {K <: Integer} = build_edges( graphs, K )

# Is the EdgeType a connecting edge between two nodes?
# Check a donor splice site by:   is_edge( edge, true )
# and an acceptor splice site by: is_edge( edge, false )
# Return: Bool
function is_edge( edge::EdgeType, left::Bool )
   eint = convert(UInt8, edge)
   eint == 0x04 && return true
   left || (eint ⊻= 0b011)
   (eint == 0x02 || eint == 0x05) && return true
   false
end

Base.intersect( arrA::Vector{SGNode}, arrB::Vector{SGNode} ) = intersect_sorted( arrA, arrB )

# Binary search for the index range of entries with node.gene == g in a gene-sorted Vector{SGNode}
@inline function gene_node_range( arr::Vector{SGNode}, g::GeneInt )
   len = length(arr)
   # findfirst >= g
   lo = 1
   hi = len + 1
   @inbounds while lo < hi
      mid = (lo + hi) >>> 1
      if arr[mid].gene < g
         lo = mid + 1
      else
         hi = mid
      end
   end
   first_idx = lo
   if first_idx > len || (@inbounds arr[first_idx].gene != g)
      return 1:0
   end
   # findlast == g (usually 1-2 elements, linear scan is faster than second binary search for small runs)
   last_idx = first_idx
   @inbounds while last_idx < len && arr[last_idx + 1].gene == g
      last_idx += 1
   end
   return first_idx:last_idx
end

# Zero-allocation equivalent of filtering intersect_sorted(arrA, arrB) to rn.gene == g:
# returns the subrange of arrB matching intersect_sorted(arrA, arrB) for gene g
@inline function intersect_sorted_gene_range( arrA::Vector{SGNode}, arrB::Vector{SGNode}, g::GeneInt )
   rB = gene_node_range(arrB, g)
   isempty(rB) && return 1:0
   rA = gene_node_range(arrA, g)
   isempty(rA) && return 1:0
   cnt = min(length(rA), length(rB))
   return first(rB):(first(rB) + cnt - 1)
end

# Intersect two sorted arrays of SGNode by gene entry
# Return the list of intersected arrB entries
function intersect_sorted( arrA::Vector{T}, arrB::Vector{T} ) where T
   res = Vector{T}()
   i,j = 1,1
   @inbounds while i <= length(arrA) && j <= length(arrB)
      if arrA[i] > arrB[j]
         j += 1
      elseif arrA[i] < arrB[j]
         i += 1
      else
         push!( res, arrB[j] )
         i += 1
         j += 1
      end
   end
   res
end


function add_kmer_edge!( kmers::Vector{SGNodeSet},
                         seq::S, l, r, left::Bool,
                         entry::SGNode ) where S <: SGSequence
   (l <= 0 || r > length(seq)) && return(zero(UInt64))
   s = copy(seq[l:r])
   ksize = r-l+1
   ind = 0 #default
   try
      curkmer = sgkmer( s )
      ind = kmer_index(curkmer)
   catch
      abstr = String(s)
      occursin( r"S|N", abstr ) && return(zero(UInt64))
      sub = replace( abstr, r"D|R", "" )
      curl,curr = l-1,r+1
      while length(sub) < ksize && curl >= 1 && curr <= length(seq)
         sub = left ? String(copy(seq[curl:curl])) * sub :
                      sub * String(copy(seq[curr:curr]))
         sub = replace( sub, r"D|R", "" )
         curr += 1
         curl -= 1
      end
      try
        curkmer = sgkmer( sub )
        ind = kmer_index(curkmer)
      catch
      end
   end
   if ind > 0
      if !isassigned(kmers, ind)
         kmers[ind] = SGNode[entry]
      else
         v = kmers[ind]
         last_el = @inbounds v[end]
         if sortlt(last_el, entry)
            push!(v, entry)
         elseif last_el.gene != entry.gene || last_el.node != entry.node
            # Rare out-of-order insertion fallback
            push!(v, entry)
            sort!(v, lt=sortlt)
            unique!(v)
         end
      end
   end
   UInt64(max(0, ind-1))
end

function build_edges( graphs::Vector{SpliceGraph}, k::Integer )
   @assert(1 <= k <= 32, "$k is out of bounds for Edges{K} set!")
   left  = Vector{SGNodeSet}(undef, 4^k)
   right = Vector{SGNodeSet}(undef, 4^k)

   for (i,g) in enumerate(graphs)
      for (j,n) in enumerate(g.nodeoffset)

         if is_edge( g.edgetype[j], true ) #&& isvalid( g.seq, (n-k-2):(n-3) ) # left edge
            lkmer = add_kmer_edge!( left, g.seq, n-k, n-1, true,  SGNode(i,j) )
            g.edgeleft[j] = reinterpret(SGKmer{k}, lkmer)
         end
         if is_edge( g.edgetype[j], false ) #&& isvalid( g.seq, n:(n+k-1) )# right edge
            rkmer = add_kmer_edge!( right, g.seq, n,  n+k-1, false, SGNode(i,j) )
            g.edgeright[j] = reinterpret(SGKmer{k}, rkmer)
         end

      end
   end
   Edges{k}(left, right)
end
