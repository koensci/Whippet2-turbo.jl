
function make_fqparser( filename; forcegzip=false )
   fopen = open( filename, "r" )
   if isgzipped( filename ) || forcegzip
      to_open = ZlibInflateInputStream( fopen, reset_on_end=true )
   else
      to_open = BufferedInputStream( fopen )
   end
   FASTQ.Reader( to_open, fill_ambiguous=DNA_A )
end

# modified for Bio v0.2 with tryread_bool!
@inline function read_chunk!( chunk, parser )
   i = 1
   while i <= length(chunk) && !eof(parser)
      read!( parser, chunk[i] )
      i += 1
   end
   while i <= length(chunk)
      pop!(chunk) # clean up if we are at the end
   end
   parser
end

function allocate_chunk( parser; size=10000 )
  chunk = Vector{eltype(parser)}( size )
  for i in 1:length(chunk)
     chunk[i] = eltype(parser)()
  end
  chunk
end

function allocate_fastq_records( size::Int=10000 )
   chunk = Vector{FASTQRecord}( undef, size )
   for i in 1:length(chunk)
      chunk[i] = FASTQRecord()
   end
   chunk
end

function process_reads!( parser, 
                        param::AlignParam, 
                        lib::SeqLibrary, 
                        quant::GraphLibQuant,
                        multi::MultiMapping{SGAlignSingle}, 
                        mod::B;
                        bufsize=1024, 
                        sam=false, 
                        qualoffset=33 ) where B <: BiasModel
   fast_lib = FastGraphLib(lib)
   if Threads.nthreads() > 1
      return _process_reads_mt!( parser, param, fast_lib, quant, multi, mod; bufsize=bufsize, sam=sam, qualoffset=qualoffset )
   end
   return _process_reads!( parser, param, fast_lib, quant, multi, mod; bufsize=bufsize, sam=sam, qualoffset=qualoffset )
end

function _process_reads!( parser, 
                          param::AlignParam, 
                          lib::FastGraphLib, 
                          quant::GraphLibQuant,
                          multi::MultiMapping{SGAlignSingle}, 
                          mod::B;
                          bufsize=1024, 
                          sam=false, 
                          qualoffset=33 ) where B <: BiasModel

   reads  = allocate_fastq_records( bufsize )
   mean_readlen = 0.0
   total        = 0
   mapped       = 0
   if sam
      stdbuf = BufferedOutputStream( stdout )
      write_sam_header( stdbuf, lib )
   end
   while length(reads) > 0
      read_chunk!( reads, parser )
      nreads = length(reads)
      total += nreads
      @inbounds for i in 1:nreads
         fill!( reads[i], qualoffset )
         align = ungapped_align( param, lib, reads[i] )
         if !isnull( align )
            biasval = count!( mod, reads[i].sequence )
            if length( align.value ) > 1
               push!( multi, align.value, biasval, quant, lib )
               sam && write_sam( stdbuf, reads[i], align.value, lib, qualoffset=qualoffset )
            else
               count!( quant, align.value[1], biasval )
               sam && write_sam( stdbuf, reads[i], align.value[1], lib, qualoffset=qualoffset )
            end
            mapped += 1
            @fastmath mean_readlen += (length(reads[i].sequence) - mean_readlen) / mapped
         end
      end
   end # end while
   if sam
      close(stdbuf)
   end
   mapped,total,mean_readlen
end

function _process_reads_mt!( parser,
                             param::AlignParam,
                             lib::FastGraphLib,
                             quant::GraphLibQuant,
                             multi::MultiMapping{SGAlignSingle},
                             mod::B;
                             bufsize=1024,
                             sam=false,
                             qualoffset=33 ) where B <: BiasModel
   nworkers = Threads.nthreads()
   chunk_sz = clamp(max(bufsize, 512 * nworkers), 2048, 8192)
   cur_reads = allocate_fastq_records( chunk_sz )
   nxt_reads = allocate_fastq_records( chunk_sz )
   alns      = Vector{Nullable{Vector{SGAlignment}}}(undef, chunk_sz)

   mean_readlen = 0.0
   total        = 0
   mapped       = 0
   if sam
      stdbuf = BufferedOutputStream( stdout )
      write_sam_header( stdbuf, lib )
   end

   read_chunk!( cur_reads, parser )
   while !isempty(cur_reads)
      nreads = length(cur_reads)
      has_more = !eof(parser)
      if has_more
         t_io = Threads.@spawn read_chunk!( nxt_reads, parser )
      end

      block = cld(nreads, nworkers)
      @sync for w in 1:nworkers
         lo = (w - 1) * block + 1
         hi = min(w * block, nreads)
         lo > hi && break
         Threads.@spawn begin
            @inbounds for i in lo:hi
               fill!( cur_reads[i], qualoffset )
               alns[i] = ungapped_align( param, lib, cur_reads[i] )
            end
         end
      end

      total += nreads
      @inbounds for i in 1:nreads
         align = alns[i]
         if !isnull( align )
            biasval = count!( mod, cur_reads[i].sequence )
            if length( align.value ) > 1
               push!( multi, align.value, biasval, quant, lib )
               sam && write_sam( stdbuf, cur_reads[i], align.value, lib, qualoffset=qualoffset )
            else
               count!( quant, align.value[1], biasval )
               sam && write_sam( stdbuf, cur_reads[i], align.value[1], lib, qualoffset=qualoffset )
            end
            mapped += 1
            @fastmath mean_readlen += (length(cur_reads[i].sequence) - mean_readlen) / mapped
         end
      end

      if has_more
         wait(t_io)
         cur_reads, nxt_reads = nxt_reads, cur_reads
      else
         break
      end
   end

   if sam
      close(stdbuf)
   end
   mapped,total,mean_readlen
end


function process_paired_reads!( fwd_parser, 
                                rev_parser, 
                                param::AlignParam,
                                lib::SeqLibrary, 
                                quant::GraphLibQuant,
                                multi::MultiMapping{SGAlignPaired}, 
                                mod::B;
                                bufsize=1024, 
                                sam=false, 
                                qualoffset=33 ) where B <: BiasModel
   fast_lib = FastGraphLib(lib)
   if Threads.nthreads() > 1
      return _process_paired_reads_mt!( fwd_parser, rev_parser, param, fast_lib, quant, multi, mod;
                                        bufsize=bufsize, sam=sam, qualoffset=qualoffset )
   end
   return _process_paired_reads!( fwd_parser, rev_parser, param, fast_lib, quant, multi, mod;
                                  bufsize=bufsize, sam=sam, qualoffset=qualoffset )
end

function _process_paired_reads!( fwd_parser, 
                                 rev_parser, 
                                 param::AlignParam,
                                 lib::FastGraphLib, 
                                 quant::GraphLibQuant,
                                 multi::MultiMapping{SGAlignPaired}, 
                                 mod::B;
                                 bufsize=1024, 
                                 sam=false, 
                                 qualoffset=33 ) where B <: BiasModel

   fwd_reads  = allocate_fastq_records( bufsize )
   rev_reads  = allocate_fastq_records( bufsize )
   fwd_buf    = Vector{Int}()
   rev_buf    = Vector{Int}()
   sizehint!(fwd_buf, param.seed_tolerate)
   sizehint!(rev_buf, param.seed_tolerate)
   mean_readlen = 0.0
   total        = 0
   mapped       = 0
   if sam
      stdbuf = BufferedOutputStream( stdout )
      write_sam_header( stdbuf, lib )
   end
   while length(fwd_reads) > 0 && length(rev_reads) > 0
      read_chunk!( fwd_reads, fwd_parser )
      read_chunk!( rev_reads, rev_parser )
      nreads = min(length(fwd_reads), length(rev_reads))
      total += nreads
      @inbounds for i in 1:nreads
         fill!( fwd_reads[i], qualoffset )
         fill!( rev_reads[i], qualoffset )
         fwd_aln,rev_aln = ungapped_align( param, lib, fwd_reads[i], rev_reads[i],
                                           fwd_buf=fwd_buf, rev_buf=rev_buf )
         if !isnull( fwd_aln ) && !isnull( rev_aln )
            biasval = count!( mod, fwd_reads[i].sequence, rev_reads[i].sequence )
            if length( fwd_aln.value ) > 1
               push!( multi, fwd_aln.value, rev_aln.value, biasval, quant, lib )
               sam && write_sam( stdbuf, fwd_reads[i], rev_reads[i], fwd_aln.value, rev_aln.value, lib,
                                 paired=true, is_pair_rc=param.is_pair_rc, qualoffset=qualoffset )
            else
               count!( quant, fwd_aln.value[1], rev_aln.value[1], biasval )
               sam && write_sam( stdbuf, fwd_reads[i], fwd_aln.value[1], lib,
                                 paired=true, fwd_mate=true, is_pair_rc=param.is_pair_rc,
                                 qualoffset=qualoffset )
               sam && write_sam( stdbuf, rev_reads[i], rev_aln.value[1], lib,
                                 paired=true, fwd_mate=false, is_pair_rc=param.is_pair_rc,
                                 qualoffset=qualoffset )
            end
            mapped += 1
            @fastmath mean_readlen += (length(fwd_reads[i].sequence) - mean_readlen) / mapped
         end
      end
   end # end while
   if sam
      close(stdbuf)
   end
   mapped,total,mean_readlen
end

function _process_paired_reads_mt!( fwd_parser,
                                    rev_parser,
                                    param::AlignParam,
                                    lib::FastGraphLib,
                                    quant::GraphLibQuant,
                                    multi::MultiMapping{SGAlignPaired},
                                    mod::B;
                                    bufsize=1024,
                                    sam=false,
                                    qualoffset=33 ) where B <: BiasModel
   nworkers = Threads.nthreads()
   chunk_sz = clamp(max(bufsize, 512 * nworkers), 2048, 8192)
   cur_fwd  = allocate_fastq_records( chunk_sz )
   cur_rev  = allocate_fastq_records( chunk_sz )
   nxt_fwd  = allocate_fastq_records( chunk_sz )
   nxt_rev  = allocate_fastq_records( chunk_sz )
   fwd_alns = Vector{Nullable{Vector{SGAlignment}}}(undef, chunk_sz)
   rev_alns = Vector{Nullable{Vector{SGAlignment}}}(undef, chunk_sz)

   w_fwd_bufs = [sizehint!(Vector{Int}(), param.seed_tolerate) for _ in 1:nworkers]
   w_rev_bufs = [sizehint!(Vector{Int}(), param.seed_tolerate) for _ in 1:nworkers]

   mean_readlen = 0.0
   total        = 0
   mapped       = 0
   if sam
      stdbuf = BufferedOutputStream( stdout )
      write_sam_header( stdbuf, lib )
   end

   # Initial read of both mate files concurrently
   t_init = Threads.@spawn read_chunk!( cur_fwd, fwd_parser )
   read_chunk!( cur_rev, rev_parser )
   wait(t_init)

   while !isempty(cur_fwd) && !isempty(cur_rev)
      nreads = min(length(cur_fwd), length(cur_rev))
      has_more = !eof(fwd_parser) && !eof(rev_parser)
      if has_more
         t_r1 = Threads.@spawn read_chunk!( nxt_fwd, fwd_parser )
         t_r2 = Threads.@spawn read_chunk!( nxt_rev, rev_parser )
      end

      block = cld(nreads, nworkers)
      @sync for w in 1:nworkers
         lo = (w - 1) * block + 1
         hi = min(w * block, nreads)
         lo > hi && break
         Threads.@spawn begin
            fb = w_fwd_bufs[w]
            rb = w_rev_bufs[w]
            @inbounds for i in lo:hi
               fill!( cur_fwd[i], qualoffset )
               fill!( cur_rev[i], qualoffset )
               fa, ra = ungapped_align( param, lib, cur_fwd[i], cur_rev[i],
                                        fwd_buf=fb, rev_buf=rb )
               fwd_alns[i] = fa
               rev_alns[i] = ra
            end
         end
      end

      total += nreads
      @inbounds for i in 1:nreads
         fwd_aln = fwd_alns[i]
         rev_aln = rev_alns[i]
         if !isnull( fwd_aln ) && !isnull( rev_aln )
            biasval = count!( mod, cur_fwd[i].sequence, cur_rev[i].sequence )
            if length( fwd_aln.value ) > 1
               push!( multi, fwd_aln.value, rev_aln.value, biasval, quant, lib )
               sam && write_sam( stdbuf, cur_fwd[i], cur_rev[i], fwd_aln.value, rev_aln.value, lib,
                                 paired=true, is_pair_rc=param.is_pair_rc, qualoffset=qualoffset )
            else
               count!( quant, fwd_aln.value[1], rev_aln.value[1], biasval )
               sam && write_sam( stdbuf, cur_fwd[i], fwd_aln.value[1], lib,
                                 paired=true, fwd_mate=true, is_pair_rc=param.is_pair_rc,
                                 qualoffset=qualoffset )
               sam && write_sam( stdbuf, cur_rev[i], rev_aln.value[1], lib,
                                 paired=true, fwd_mate=false, is_pair_rc=param.is_pair_rc,
                                 qualoffset=qualoffset )
            end
            mapped += 1
            @fastmath mean_readlen += (length(cur_fwd[i].sequence) - mean_readlen) / mapped
         end
      end

      if has_more
         wait(t_r1)
         wait(t_r2)
         cur_fwd, nxt_fwd = nxt_fwd, cur_fwd
         cur_rev, nxt_rev = nxt_rev, cur_rev
      else
         break
      end
   end

   if sam
      close(stdbuf)
   end
   mapped,total,mean_readlen
end
