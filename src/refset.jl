
# Single RefTx entry
struct RefTx
   info::TxInfo
   don::CoordTuple
   acc::CoordTuple
   length::Float64
end

# Single RefGene entry
struct RefGene
   info::GeneInfo
   don::CoordTuple
   acc::CoordTuple
   txst::CoordTuple
   txen::CoordTuple
   exons::CoordTree
   noveldon::CoordTuple
   novelacc::CoordTuple
   exonexpr::Float64
   length::Float64
   reftx::Vector{RefTx}
end

const RefSet = Dict{GeneName,RefGene}

# this guy maps an array of strings or substrings to Int parsing and then into a tuple
# which is then sliced based on the l (left) and r (right) adjusters, c is the mathmatical adjuster to the data
parse_splice( subarray; l=0, r=0, c=0 ) = tuple(map( x->convert(CoordInt,parse(Int,x)+c), subarray )...)[(1+l):(end-r)]
parse_coordint( i ) = convert(CoordInt, parse(Int, i))
parse_coordtup( i; c=0 ) = tuple(parse_coordint(i)+CoordInt(c))

unique_tuple( tup1::Tuple{}, tup2::Tuple{} ) = ()
unique_tuple( tup1::Tuple{Vararg{T}}, tup2::Tuple{} ) where T = unique_tuple(tup1, tup1)
unique_tuple( tup1::Tuple{}, tup2::Tuple{Vararg{T}} ) where T = unique_tuple(tup2, tup2)

function unique_tuple( tup1::Tuple{Vararg{T}}, tup2::Tuple{Vararg{T}}) where T
   uniq = SortedSet{T,Base.Order.ForwardOrdering}()
   for i in tup1
      push!(uniq, i)
   end
   for j in tup2
      push!(uniq, j)
   end
   tuple(uniq...)
end

function minimum_threshold!( dict::Dict{K,V}, limit::V, range::UnitRange ) where {K, V <: Number}
   for k in collect(keys(dict))
      if dict[k] < limit || !(k in range)
         delete!( dict, k )
      end
   end
end

function load_refflat( fh; txbool=true )
   error("ERROR: Annotation files in refflat format have been deprecated as of Whippet v0.11, please use --gtf!")
end

function fetch_meta( var::String, meta; off=1 )
   val = "_" # empty
   for i in off:2:length(meta)
      if meta[i] == var && i <= length(meta)
         return string(meta[i+1]),i+1
      end
   end
   val,0
end

# Fast single-pass GTF attribute extractor (zero intermediate array allocations)
function _extract_gtf_attrs( meta::AbstractString )
   geneid  = "_"
   tranid  = "_"
   support = "_"
   tag     = "_"
   n = lastindex(meta)
   i = firstindex(meta)
   @inbounds while i <= n
      # skip spaces/semicolons/tabs
      while i <= n
         c = meta[i]
         (c == ' ' || c == ';' || c == '\t') || break
         i = nextind(meta, i)
      end
      i > n && break
      # read key
      kstart = i
      while i <= n
         c = meta[i]
         (c == ' ' || c == '\t' || c == ';') && break
         i = nextind(meta, i)
      end
      kend = prevind(meta, i)
      # skip spaces before value
      while i <= n && (meta[i] == ' ' || meta[i] == '\t')
         i = nextind(meta, i)
      end
      i > n && break
      if meta[i] == ';'
         i = nextind(meta, i)
         continue
      end
      # read value (quoted or unquoted)
      if meta[i] == '"'
         vstart = nextind(meta, i)
         i = vstart
         while i <= n && meta[i] != '"'
            i = nextind(meta, i)
         end
         vend = prevind(meta, i)
         i <= n && (i = nextind(meta, i))
      else
         vstart = i
         while i <= n && meta[i] != ';' && meta[i] != ' ' && meta[i] != '\t'
            i = nextind(meta, i)
         end
         vend = prevind(meta, i)
      end
      key = SubString(meta, kstart, kend)
      if key == "gene_id" && geneid == "_"
         geneid = vstart <= vend ? String(SubString(meta, vstart, vend)) : ""
      elseif key == "transcript_id" && tranid == "_"
         tranid = vstart <= vend ? String(SubString(meta, vstart, vend)) : ""
      elseif key == "transcript_support_level" && support == "_"
         support = vstart <= vend ? String(SubString(meta, vstart, vend)) : ""
      elseif key == "tag" && tag == "_"
         tag = vstart <= vend ? String(SubString(meta, vstart, vend)) : ""
      end
   end
   return geneid, tranid, support, tag
end

function load_gtf( fh; txbool=true, 
                       suppress=false,
                       basic_or_ssonly=false,
                       usebam=false, 
                       bamreader=Nullable{BAM.Reader}(), 
                       bamreads=2, 
                       bamoneknown=false )

   # Accumulate coordinates in flat Vectors per gene (avoids O(N*M) tuple specialization)
   gninfo   = Dict{GeneName,GeneInfo}()
   gndon    = Dict{GeneName,Vector{CoordInt}}()
   gnacc    = Dict{GeneName,Vector{CoordInt}}()
   gntxst   = Dict{GeneName,Vector{CoordInt}}()
   gntxen   = Dict{GeneName,Vector{CoordInt}}()
   gnexons  = Dict{GeneName,CoordTree}()
   gnreftx  = Dict{GeneName,Vector{RefTx}}()
   gnlen    = Dict{GeneName,Float64}()
   gncnt    = Dict{GeneName,Int}()

   # RefSet variables
   geneset  = Dict{GeneName,RefGene}()

   used_txn = Dict{GeneName,Bool}()
   warning  = false
   warning_num = 0
   ssonly   = false

   txnum = 75000
   sizehint!(geneset, txnum >> 1)
   sizehint!(gninfo, txnum >> 1)
   sizehint!(gndon, txnum >> 1)
   sizehint!(gnacc, txnum >> 1)
   sizehint!(gntxst, txnum >> 1)
   sizehint!(gntxen, txnum >> 1)
   sizehint!(gnexons, txnum >> 1)
   sizehint!(gnreftx, txnum >> 1)

   curtran = ""
   curgene = ""
   curchrom = ""
   curstran = '+'
   curbasic = false

   trandon = Vector{CoordInt}()
   tranacc = Vector{CoordInt}()
   txlen   = 0

   function private_add_transcript!( curtran, curgene, curchrom, curstran, trandon, tranacc, txlen;
                                     splice_sites_only=false )
      isempty(tranacc) && return
      sort!(trandon)
      sort!(tranacc)

      txS = first(tranacc)
      txE = last(trandon)

      txinfo = TxInfo(curtran, txS, txE, UInt16(length(tranacc)))

      has_info = haskey(gninfo, curgene)
      if has_info
         curchrom == gninfo[curgene].name || return # can have only one chrom
      end

      dvec = get!(Vector{CoordInt}, gndon, curgene)
      avec = get!(Vector{CoordInt}, gnacc, curgene)
      nt = length(trandon)
      if nt > 1
         @inbounds for k in 1:(nt - 1)
            push!(dvec, trandon[k])
         end
         @inbounds for k in 2:nt
            push!(avec, tranacc[k])
         end
      end

      if !splice_sites_only
         if has_info
            push!(gntxst[curgene], txS)
            push!(gntxen[curgene], txE)
            gnlen[curgene] += txlen
            gncnt[curgene] += 1
            push!(gnreftx[curgene], RefTx(txinfo, Tuple(trandon), Tuple(tranacc), txlen))
         else
            gninfo[curgene]  = GeneInfo(curgene, string(curchrom), curstran)
            gntxst[curgene]  = CoordInt[txS]
            gntxen[curgene]  = CoordInt[txE]
            gnlen[curgene]   = Float64(txlen)
            gncnt[curgene]   = 1
            gnreftx[curgene] = RefTx[RefTx(txinfo, Tuple(trandon), Tuple(tranacc), txlen)]
         end
      end
   end

   for raw_l in eachline(fh)
      isempty(raw_l) && continue
      l = (!isempty(raw_l) && raw_l[end] == '\r') ? SubString(raw_l, 1, prevind(raw_l, lastindex(raw_l))) : raw_l
      (isempty(l) || l[1] == '#') && continue # ignore comment lines

      t1 = findfirst('\t', l); t1 === nothing && continue
      t2 = findnext('\t', l, t1 + 1); t2 === nothing && continue
      t3 = findnext('\t', l, t2 + 1); t3 === nothing && continue
      # Fast filter on column 3 ("exon") before parsing remaining columns or attributes
      SubString(l, t2 + 1, t3 - 1) == "exon" || continue

      t4 = findnext('\t', l, t3 + 1); t4 === nothing && continue
      t5 = findnext('\t', l, t4 + 1); t5 === nothing && continue
      t6 = findnext('\t', l, t5 + 1); t6 === nothing && continue
      t7 = findnext('\t', l, t6 + 1); t7 === nothing && continue
      t8 = findnext('\t', l, t7 + 1); t8 === nothing && continue

      chrom  = SubString(l, 1, t1 - 1)
      st     = SubString(l, t3 + 1, t4 - 1)
      en     = SubString(l, t4 + 1, t5 - 1)
      strand = SubString(l, t6 + 1, t7 - 1)
      meta   = SubString(l, t8 + 1)

      geneid, tranid, support, tag = _extract_gtf_attrs(meta)

      # if we observe low transcript support levels, we will warn the user
      # unless we are using gencode system of basic tx for txS and txE and others
      # only for splice site support
      if !basic_or_ssonly && 
         (support != "_" && support != "1" && support != "2" &&  support != "NA")

         if suppress
            continue
         elseif !warning
            println(stderr, "")
            @warn("Using low quality Transcript Support Levels (TSL 3+) in your GTF file is not recommended!\nFor more information on TSL, see: http://www.ensembl.org/Help/Glossary?id=492\n\nIf you would like Whippet to ignore these when building its index, use `--suppress-low-tsl` option!\n\n")

            warning=true
         end

      elseif haskey(used_txn, tranid)
         error("ERROR: GTF file is not in valid GTF2.2 format!\n\nERROR: Annotation entries for 'transcript_id' $tranid has already been fully processed and closed.\nHINT: All GTF lines with the same 'transcript_id' must be adjacent in the GTF file and referring to the same transcript and gene!")
      elseif tranid == geneid && tranid != curtran
         if warning_num < 25
            @warn("Generally 'transcript_id' should not equal 'gene_id' but does at $tranid == $geneid;")
         elseif warning_num == 25
            @warn("... NOTE: 'transcript_id' == 'gene_id' will work OK for single isoform genes, but will not produce expected behavior for multi-isoform genes!\n... similar warnings will be suppressed; disregard if 'transcript_id' == 'gene_id' is intentional\n")
         end
         warning_num += 1
      end

      if tranid != curtran
         # add transcript and clean up
         ssonly = basic_or_ssonly && !curbasic ? true : false
         curtran != "" && private_add_transcript!( curtran, curgene, curchrom, curstran, trandon, tranacc, txlen, splice_sites_only=ssonly )

         used_txn[curtran] = true

         curtran = tranid
         curgene = geneid
         curchrom = String(chrom)
         curstran = strand[1]
         curbasic = tag == "basic" ? true : false

         empty!(trandon)
         empty!(tranacc)
         txlen = 0
      end

      sval = parse_coordint(st)
      eval = parse_coordint(en)

      push!(tranacc, sval)
      push!(trandon, eval)
      txlen += eval - sval + 1

      insval = IntervalTrees.Interval{CoordInt}(sval, eval)
      ex_tree = get!(CoordTree, gnexons, curgene)
      if !haskey(ex_tree, (sval, eval))
         push!(ex_tree, insval)
      end
   end

   ssonly = basic_or_ssonly && !curbasic ? true : false
   curtran != "" && private_add_transcript!( curtran, curgene, curchrom, curstran, trandon, tranacc, txlen, splice_sites_only=ssonly )

   bamnames = usebam ? Set(map(x->values(x)[1], findall(BAM.header(get(bamreader)), "SQ"))) : Set{String}()

   novel_ss = 0
   annot_ss = 0
   # iterate through genes, if usebam, add novel splice sites and calculate relative exonexpr
   # set RefGene at the end of each iteration
   for gene in keys(gninfo)

      novelacc = Dict{CoordInt,Int}()
      noveldon = Dict{CoordInt,Int}()
      exonexpr = 0.0
      meanleng = gnlen[gene] / gncnt[gene]
      seqname  = gninfo[gene].name
      strand   = gninfo[gene].strand
      bamwasused = false

      don_vec  = gndon[gene];  sort!(don_vec);  unique!(don_vec)
      acc_vec  = gnacc[gene];  sort!(acc_vec);  unique!(acc_vec)
      txst_vec = gntxst[gene]; sort!(txst_vec); unique!(txst_vec)
      txen_vec = gntxen[gene]; sort!(txen_vec); unique!(txen_vec)

      generange = Int64(first(txst_vec)):Int64(last(txen_vec))

      annoacc_vec = copy(acc_vec)
      annodon_vec = copy(don_vec)
      novelacc_vec = CoordInt[]
      noveldon_vec = CoordInt[]

      if usebam && seqname in bamnames
         bamwasused = true
         exoncount = process_records!( get(bamreader), seqname, generange, strand,
                                       gnexons[gene], union(acc_vec, don_vec), bamoneknown,
                                       novelacc, noveldon )
         exonexpr  = exoncount / (meanleng - 100)

         # clean up cryptic splice-sites
         minimum_threshold!( novelacc, bamreads, generange )
         minimum_threshold!( noveldon, bamreads, generange )
         novelacc_vec = CoordInt[CoordInt(k) for k in keys(novelacc)]
         noveldon_vec = CoordInt[CoordInt(k) for k in keys(noveldon)]
         sort!(novelacc_vec); sort!(noveldon_vec)

         # add splice-sites from bam
         n_before = length(acc_vec) + length(don_vec)
         append!(acc_vec, novelacc_vec); sort!(acc_vec); unique!(acc_vec)
         append!(don_vec, noveldon_vec); sort!(don_vec); unique!(don_vec)
         novel_ss += (length(acc_vec) + length(don_vec)) - n_before
      end

      # clean up any txst or txen that match a splice site
      txst_clean = setdiff(txst_vec, acc_vec)
      txen_clean = setdiff(txen_vec, don_vec)

      don_clean  = filter!(x -> x in generange, don_vec)
      acc_clean  = filter!(x -> x in generange, acc_vec)

      annot_ss += length(don_clean) + length(acc_clean)

      geneset[gene] = RefGene( gninfo[gene],
                               Tuple(don_clean),
                               Tuple(acc_clean),
                               Tuple(txst_clean),
                               Tuple(txen_clean),
                               gnexons[gene],
                               bamwasused ? Tuple(setdiff(noveldon_vec, annodon_vec)) : Tuple{}(),
                               bamwasused ? Tuple(setdiff(novelacc_vec, annoacc_vec)) : Tuple{}(),
                               exonexpr,
                               meanleng,
                               gnreftx[gene] )
   end

   if usebam
      println(stderr, "Found $novel_ss new splice-sites from BAM file..")
   end
   println(stderr, "Loaded $annot_ss annotated splice-sites from GTF file..")

   return geneset
end
