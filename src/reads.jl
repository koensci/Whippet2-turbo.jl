
const _FAST_GZIP_PROG = Ref{Union{Nothing, String}}(nothing)
const _FAST_GZIP_CHECKED = Ref{Bool}(false)

function _fast_gzip_prog()
   if !_FAST_GZIP_CHECKED[]
      p = nothing
      if Sys.isunix()
         p = Sys.which("igzip")
         if p === nothing
            p = Sys.which("pigz")
         end
      end
      _FAST_GZIP_PROG[] = p
      _FAST_GZIP_CHECKED[] = true
   end
   _FAST_GZIP_PROG[]
end

# Direct POSIX popen/fread IO wrapper to avoid libuv iolock contention on multi-threaded reads.
mutable struct _PopenReadStream <: IO
   fp::Ptr{Cvoid}
   closed::Bool
   function _PopenReadStream(cmd::String)
      fp = ccall(:popen, Ptr{Cvoid}, (Cstring, Cstring), cmd, "r")
      fp == C_NULL && error("popen failed for command: $cmd")
      s = new(fp, false)
      finalizer(close, s)
      return s
   end
end

Base.isopen(s::_PopenReadStream) = !s.closed
Base.eof(s::_PopenReadStream) = s.closed || (ccall(:feof, Cint, (Ptr{Cvoid},), s.fp) != 0)

function Base.readbytes!(s::_PopenReadStream, b::AbstractVector{UInt8}, nb=length(b))
   s.closed && return 0
   n = min(Int(nb), length(b))
   n <= 0 && return 0
   got = GC.@preserve b ccall(:fread, Csize_t, (Ptr{UInt8}, Csize_t, Csize_t, Ptr{Cvoid}), pointer(b), 1, n, s.fp)
   return Int(got)
end

function BufferedStreams.readbytes!(s::_PopenReadStream, b::Vector{UInt8}, from::Int, to::Int)
   s.closed && return 0
   n = to - from + 1
   n <= 0 && return 0
   got = GC.@preserve b ccall(:fread, Csize_t, (Ptr{UInt8}, Csize_t, Csize_t, Ptr{Cvoid}), pointer(b, from), 1, n, s.fp)
   return Int(got)
end

function Base.close(s::_PopenReadStream)
   if !s.closed && s.fp != C_NULL
      ccall(:pclose, Cint, (Ptr{Cvoid},), s.fp)
      s.fp = C_NULL
      s.closed = true
   end
   nothing
end

@inline _sh_quote(s::AbstractString) = "'" * replace(String(s), "'" => "'\\''") * "'"

function make_fqparser( filename; forcegzip=false )
   if isgzipped( filename ) || forcegzip
      prog = isfile( filename ) ? _fast_gzip_prog() : nothing
      if prog !== nothing
         cmd = _sh_quote(prog) * " -dc " * _sh_quote(filename)
         to_open = BufferedInputStream( _PopenReadStream(cmd), 262144 )
      else
         fopen = open( filename, "r" )
         to_open = ZlibInflateInputStream( fopen, reset_on_end=true, bufsize=262144 )
      end
   else
      fopen = open( filename, "r" )
      to_open = BufferedInputStream( fopen, 262144 )
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

@inline function read_chunk_n!( chunk::Vector{FASTQRecord}, parser )::Int
   i = 0
   cap = length(chunk)
   @inbounds while i < cap && !eof(parser)
      i += 1
      read!( parser, chunk[i] )
   end
   i
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

struct PairedChunkBuf
   fwd::Vector{FASTQRecord}
   rev::Vector{FASTQRecord}
   fwd_aln::Vector{SGAlignment}
   rev_aln::Vector{SGAlignment}
   fwd_multi::Vector{Vector{SGAlignment}}
   rev_multi::Vector{Vector{SGAlignment}}
   nhits::Vector{Int32}
end

function PairedChunkBuf( sz::Int )
   fwd = allocate_fastq_records( sz )
   rev = allocate_fastq_records( sz )
   fwd_aln = [SGAlignment(zero(UInt32), zero(ReadLengthInt), sizehint!(SGAlignNode[], 4), true, false) for _ in 1:sz]
   rev_aln = [SGAlignment(zero(UInt32), zero(ReadLengthInt), sizehint!(SGAlignNode[], 4), true, false) for _ in 1:sz]
   fwd_multi = [SGAlignment[] for _ in 1:sz]
   rev_multi = [SGAlignment[] for _ in 1:sz]
   nhits = zeros(Int32, sz)
   PairedChunkBuf( fwd, rev, fwd_aln, rev_aln, fwd_multi, rev_multi, nhits )
end

function _align_paired_chunk_mt!( buf::PairedChunkBuf,
                                  nreads::Int,
                                  param::AlignParam,
                                  lib::FastGraphLib,
                                  qualoffset::Int,
                                  nworkers::Int,
                                  w_fwd_bufs::Vector{Vector{Int}},
                                  w_rev_bufs::Vector{Vector{Int}} )
   block = cld(nreads, nworkers)
   @sync for w in 1:nworkers
      lo = (w - 1) * block + 1
      hi = min(w * block, nreads)
      lo > hi && break
      Threads.@spawn begin
         fb = @inbounds w_fwd_bufs[w]
         rb = @inbounds w_rev_bufs[w]
         fwd = buf.fwd
         rev = buf.rev
         fwd_aln = buf.fwd_aln
         rev_aln = buf.rev_aln
         fwd_multi = buf.fwd_multi
         rev_multi = buf.rev_multi
         nhits = buf.nhits
         @inbounds for i in lo:hi
            fill!( fwd[i], qualoffset )
            fill!( rev[i], qualoffset )
            nhits[i] = ungapped_align_into!( param, lib, fwd[i], rev[i],
                                             fwd_aln, rev_aln,
                                             fwd_multi[i], rev_multi[i], i;
                                             fwd_buf=fb, rev_buf=rb )
         end
      end
   end
end

function _count_paired_chunk!( buf::PairedChunkBuf,
                               nreads::Int,
                               param::AlignParam,
                               lib::FastGraphLib,
                               quant::GraphLibQuant,
                               multi::MultiMapping{SGAlignPaired},
                               mod::B,
                               sam::Bool,
                               stdbuf,
                               qualoffset::Int,
                               mapped::Int,
                               total::Int,
                               mean_readlen::Float64 ) where B <: BiasModel
   fwd = buf.fwd
   rev = buf.rev
   fwd_aln = buf.fwd_aln
   rev_aln = buf.rev_aln
   fwd_multi = buf.fwd_multi
   rev_multi = buf.rev_multi
   nhits = buf.nhits
   total += nreads
   @inbounds for i in 1:nreads
      nh = nhits[i]
      if nh > Int32(0)
         biasval = count!( mod, fwd[i].sequence, rev[i].sequence )
         if nh > Int32(1)
            push!( multi, fwd_multi[i], rev_multi[i], biasval, quant, lib )
            sam && write_sam( stdbuf, fwd[i], rev[i], fwd_multi[i], rev_multi[i], lib,
                              paired=true, is_pair_rc=param.is_pair_rc, qualoffset=qualoffset )
         else
            count!( quant, fwd_aln[i], rev_aln[i], biasval )
            sam && write_sam( stdbuf, fwd[i], fwd_aln[i], lib,
                              paired=true, fwd_mate=true, is_pair_rc=param.is_pair_rc,
                              qualoffset=qualoffset )
            sam && write_sam( stdbuf, rev[i], rev_aln[i], lib,
                              paired=true, fwd_mate=false, is_pair_rc=param.is_pair_rc,
                              qualoffset=qualoffset )
         end
         mapped += 1
         @fastmath mean_readlen += (length(fwd[i].sequence) - mean_readlen) / mapped
      end
   end
   return mapped, total, mean_readlen
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
   fwd_buf    = sizehint!(Vector{Int}(), param.seed_tolerate)
   rev_buf    = sizehint!(Vector{Int}(), param.seed_tolerate)
   fwd_alns   = [SGAlignment(zero(UInt32), zero(ReadLengthInt), sizehint!(SGAlignNode[], 4), true, false)]
   rev_alns   = [SGAlignment(zero(UInt32), zero(ReadLengthInt), sizehint!(SGAlignNode[], 4), true, false)]
   fwd_multi  = SGAlignment[]
   rev_multi  = SGAlignment[]
   mean_readlen = 0.0
   total        = 0
   mapped       = 0
   stdbuf = sam ? BufferedOutputStream( stdout ) : nothing
   sam && write_sam_header( stdbuf, lib )
   while true
      n1 = read_chunk_n!( fwd_reads, fwd_parser )
      n2 = read_chunk_n!( rev_reads, rev_parser )
      nreads = min(n1, n2)
      nreads == 0 && break
      total += nreads
      @inbounds for i in 1:nreads
         fill!( fwd_reads[i], qualoffset )
         fill!( rev_reads[i], qualoffset )
         nh = ungapped_align_into!( param, lib, fwd_reads[i], rev_reads[i],
                                    fwd_alns, rev_alns, fwd_multi, rev_multi, 1;
                                    fwd_buf=fwd_buf, rev_buf=rev_buf )
         if nh > Int32(0)
            biasval = count!( mod, fwd_reads[i].sequence, rev_reads[i].sequence )
            if nh > Int32(1)
               push!( multi, fwd_multi, rev_multi, biasval, quant, lib )
               sam && write_sam( stdbuf, fwd_reads[i], rev_reads[i], fwd_multi, rev_multi, lib,
                                 paired=true, is_pair_rc=param.is_pair_rc, qualoffset=qualoffset )
            else
               count!( quant, fwd_alns[1], rev_alns[1], biasval )
               sam && write_sam( stdbuf, fwd_reads[i], fwd_alns[1], lib,
                                 paired=true, fwd_mate=true, is_pair_rc=param.is_pair_rc,
                                 qualoffset=qualoffset )
               sam && write_sam( stdbuf, rev_reads[i], rev_alns[1], lib,
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
   chunk_sz = clamp(max(bufsize, 256 * nworkers), 2048, 4096)

   # Three disjoint chunk buffers for 3-stage pipeline:
   #   Stage 1 (b_io):  concurrent R1/R2 FASTQ reading
   #   Stage 2 (b_aln): parallel alignment across workers
   #   Stage 3 (b_cnt): serial deterministic counting, bias update, and optional SAM output
   b_aln = PairedChunkBuf( chunk_sz )
   b_io  = PairedChunkBuf( chunk_sz )
   b_cnt = PairedChunkBuf( chunk_sz )

   w_fwd_bufs = [sizehint!(Vector{Int}(), param.seed_tolerate) for _ in 1:nworkers]
   w_rev_bufs = [sizehint!(Vector{Int}(), param.seed_tolerate) for _ in 1:nworkers]

   mean_readlen = 0.0
   total        = 0
   mapped       = 0
   stdbuf = sam ? BufferedOutputStream( stdout ) : nothing
   sam && write_sam_header( stdbuf, lib )

   # Bootstrap: read first chunk concurrently into b_aln
   t_init = Threads.@spawn read_chunk_n!( b_aln.fwd, fwd_parser )
   n_r2   = read_chunk_n!( b_aln.rev, rev_parser )
   n_aln  = min(fetch(t_init)::Int, n_r2)
   n_cnt  = 0

   while n_aln > 0 || n_cnt > 0
      do_io = (n_aln > 0) && !eof(fwd_parser) && !eof(rev_parser)
      if do_io
         t_io1 = Threads.@spawn read_chunk_n!( b_io.fwd, fwd_parser )
         t_io2 = Threads.@spawn read_chunk_n!( b_io.rev, rev_parser )
      end

      if n_cnt > 0
         t_cnt = Threads.@spawn _count_paired_chunk!( b_cnt, n_cnt, param, lib, quant, multi, mod,
                                                      sam, stdbuf, qualoffset,
                                                      mapped, total, mean_readlen )
      end

      if n_aln > 0
         _align_paired_chunk_mt!( b_aln, n_aln, param, lib, qualoffset,
                                  nworkers, w_fwd_bufs, w_rev_bufs )
      end

      if n_cnt > 0
         mapped, total, mean_readlen = fetch(t_cnt)::Tuple{Int, Int, Float64}
      end

      n_next = 0
      if do_io
         n_next = min(fetch(t_io1)::Int, fetch(t_io2)::Int)
      end

      b_cnt, b_aln, b_io = b_aln, b_io, b_cnt
      n_cnt = n_aln
      n_aln = n_next
   end

   if sam
      close(stdbuf)
   end
   mapped,total,mean_readlen
end
