
abstract type ReadCounter end
abstract type BiasModel end

function pushto!(t::IntervalTrees.IntervalBTree{K, V, B}, key::I, value) where {K, V, B, I <: AbstractInterval{K}}
    return _push!(t.root, key, value)
end

function _push!(t::IntervalTrees.InternalNode{K, V, B}, key::AbstractInterval{K}, value) where {K, V, B}
    i = IntervalTrees.findidx(t, key)
    if 1 <= length(t) - 1 && key >= t.keys[i]
        return _push!(t.children[i+1], key, value)
    else
        return _push!(t.children[i], key, value)
    end
end

function _push!(t::IntervalTrees.LeafNode{K, V, B}, key::AbstractInterval{K}, value) where {K, V, B}
    i = IntervalTrees.findidx(t, key)
    if 1 <= i <= length(t) &&
        first(t.entries[i]) == first(key) && last(t.entries[i]) == last(key)
        push!( t.entries[i].value, value )
        return true
    end
    return false
end

function pushzero!( collection::Dict{K,V}, key::K, value ) where {K, V}
   if !haskey( collection, key )
      collection[key] = zero(V)
   end
   push!( collection[key], value )
end

function pushzero!( collection::IntervalMap{K,V}, key, value ) where {K <: NodeNum, V}
   if !haskey( collection, (key.first,key.last) )
      collection[key] = zero(V)
   end
   pushto!( collection, key, value )
end

function Base.get( cnt::R ) where R <: ReadCounter
   if !cnt.isadjusted
      error("Cannot get() read count from model thats not already adjusted!")
   else
      return cnt.count
   end
end

function Base.push!( cnt::R, value::F ) where {R <: ReadCounter, F <: AbstractFloat}
   if cnt.isadjusted
      cnt.count += value
   else
      cnt.count += value
      cnt.isadjusted = true
      #error("Cannot push! a float value to PrimerBiasCounter unless isadjusted is true!")
   end
   cnt
end

Base.identity( cnt::R ) where R <: ReadCounter = get( cnt )
Base.identity( cnt::IntervalValue{I,R} ) where R <: ReadCounter = get( cnt.value )

struct DefaultBiasMod <: BiasModel end

mutable struct DefaultCounter <: ReadCounter
   count::Float64
   isadjusted::Bool

   DefaultCounter() = new( 1.0, true )
   DefaultCounter( i::Float64 ) = new( i, true )
end

const DEFAULTCOUNTER_ZERO = DefaultCounter(0.0)
const DEFAULTCOUNTER_ONE  = DefaultCounter(1.0)

Base.zero( ::Type{DefaultCounter} ) = DefaultCounter(0.0) #DEFAULTCOUNTER_ZERO
Base.one( ::Type{DefaultCounter} )  = DefaultCounter(1.0) #DEFAULTCOUNTER_ONE
Base.push!( cnt::DefaultCounter, value::Float64 ) = (cnt.count += value)
Base.push!( cnt::DefaultCounter, value::DefaultCounter ) = (cnt.count += get(value))

primer_adjust!( cnt::DefaultCounter, m::B ) where B <: BiasModel = (cnt.isadjusted = true)
gc_adjust!( cnt::DefaultCounter, m::B ) where B <: BiasModel = (cnt.isadjusted = true)

count!( mod::DefaultBiasMod, _ ) = DEFAULTCOUNTER_ONE
count!( mod::DefaultBiasMod, x, y ) = DEFAULTCOUNTER_ONE

# 5' Sequence bias model
# Implementation of Hansen et al, NAR 2010 

# hash by identity for speed
Base.hash( v::UInt16 ) = convert(UInt64, v)

function add_to_back!( back::Vector{Float64}, exp::ExpectedGC, weight::Float64 )
   @inbounds for i in 1:length(exp)
      back[i] += exp[i] * weight
   end
   back
end

primer_normalize!( mod::DefaultBiasMod )   = 1.0
gc_normalize!( mod::DefaultBiasMod, x, y ) = 1.0

struct JointBiasTemp
   primer::Vector{UInt16}
   gc::Vector{Int32}
end

struct JointBiasMod <: BiasModel
   fore::Vector{Float64}
   back::Vector{Float64}
   size::Int64
   backoffset::Int64
   backnpos::Int64
   foreoffset::Int64
   forenpos::Int64

   gcfore::Vector{Float64}
   gcback::Vector{Float64}
   gcsize::Int64
   gcwidth::Float64
   gcoffset::Int64
   gcincrement::Int64

   temp::JointBiasTemp

   function JointBiasMod( ksize::Int=6, 
                          backoff::Int=23, 
                          backn::Int=10, 
                          foreoff::Int=1, 
                          foren::Int=2,
                          gcoffset::Int=5, 
                          gcincrement::Int=10; 
                          gc_param = GCBiasParams(50, 0.05) )

      fore = zeros( 4^ksize )
      back = zeros( 4^ksize )
      primer_temp = zeros(UInt16, foren)

      gcfore = zeros(Float64, Int(ceil(1.0 / gc_param.width)+1))
      gcback = zeros(Float64, Int(ceil(1.0 / gc_param.width)+1))
      gctemp = zeros(Int32, 100)

      temp = JointBiasTemp( primer_temp, gctemp )

      new( fore, back, ksize, backoff, backn, foreoff, foren,
           gcfore, gcback, gc_param.length, gc_param.width, gcoffset, gcincrement, 
           temp )
   end
end

function primer_normalize!( mod::JointBiasMod )
   foresum = sum(mod.fore)
   backsum = sum(mod.back)
   @inbounds for i in 1:length(mod.fore)
      if foresum > 0.0
         @fastmath mod.fore[i] = mod.fore[i] / foresum
      end
      if backsum > 0.0
         @fastmath mod.back[i] = mod.back[i] / backsum
      end
   end
end

function primer_count!( mod::JointBiasMod, seq::LongSequence{DNAAlphabet{2}} )
   k = mod.size
   seqlen = length(seq)
   mask = UInt16((1 << (2 * k)) - 1)
   fore = mod.fore
   back = mod.back
   temp_primer = mod.temp.primer

   # 5 prime foreground bias (rolling 2-bit k-mer)
   foreoff = mod.foreoffset
   foren = min(mod.forenpos, max(0, seqlen - foreoff - k + 2))
   if foren > 0
      hept = kmer_index_trailing_sub(UInt16, seq, foreoff, foreoff + k - 1)
      @inbounds begin
         fore[Int(hept) + 1] += 1.0
         temp_primer[1] = hept
         for idx in 2:foren
            hept = ((hept << 2) | UInt16(raw_getindex(seq, foreoff + idx + k - 2))) & mask
            fore[Int(hept) + 1] += 1.0
            temp_primer[idx] = hept
         end
      end
   end

   # background bias (rolling 2-bit k-mer)
   backoff = mod.backoffset
   backn = min(mod.backnpos, max(0, seqlen - backoff - k + 2))
   if backn > 0
      hept = kmer_index_trailing_sub(UInt16, seq, backoff, backoff + k - 1)
      @inbounds begin
         back[Int(hept) + 1] += 1.0
         for idx in 2:backn
            hept = ((hept << 2) | UInt16(raw_getindex(seq, backoff + idx + k - 2))) & mask
            back[Int(hept) + 1] += 1.0
         end
      end
   end

   return mod.temp
end

function primer_count!( mod::JointBiasMod, seq::BioSequence{A} ) where A <: BioSequences.Alphabet
   # 5 prime bias
   idx = 1
   k = mod.size
   seqlen = length(seq)
   fore = mod.fore
   back = mod.back
   temp_primer = mod.temp.primer
   fore_stop = min(mod.foreoffset + mod.forenpos - 1, seqlen - k + 1)
   @inbounds for i in mod.foreoffset:fore_stop
      hept = kmer_index_trailing_sub(UInt16, seq, i, i + k - 1)
      fore[Int(hept) + 1] += 1.0
      temp_primer[idx] = hept
      idx += 1
   end
   back_stop = min(mod.backoffset + mod.backnpos - 1, seqlen - k + 1)
   @inbounds for i in mod.backoffset:back_stop
      back[Int(kmer_index_trailing_sub(UInt16, seq, i, i + k - 1)) + 1] += 1.0
   end

   return mod.temp
end

@inline _is_gc_2bit( b::UInt8 ) = (b == 0x01) | (b == 0x02)

@inline function _push_gc_temp!( temp_gc::Vector{Int32}, idx::Int, bin::Int )
   if idx > length(temp_gc)
      oldlen = length(temp_gc)
      resize!(temp_gc, max(idx, oldlen * 2))
      @inbounds for k in (oldlen + 1):length(temp_gc)
         temp_gc[k] = zero(Int32)
      end
   end
   @inbounds temp_gc[idx] = Int32(bin)
   return idx + 1
end

function _gc_count!( mod::JointBiasMod, seq::LongSequence{DNAAlphabet{2}}, start_idx::Int )
   # GC fragment bias (zero-allocation sliding window)
   gcoff = mod.gcoffset
   gcinc = mod.gcincrement
   gcsize = mod.gcsize
   gcwidth = mod.gcwidth
   gcfore = mod.gcfore
   temp_gc = mod.temp.gc
   stop_i = length(seq) - gcsize
   gcoff <= stop_i || return start_idx

   inv_gcsize = Float64(gcsize)
   idx = start_idx
   if gcinc < gcsize
      gc_cnt = 0
      @inbounds for pos in gcoff:(gcoff + gcsize - 1)
         gc_cnt += ifelse(_is_gc_2bit(raw_getindex(seq, pos)), 1, 0)
      end
      bin = Int(div(gc_cnt / inv_gcsize, gcwidth)) + 1
      @inbounds gcfore[bin] += 1.0
      idx = _push_gc_temp!(temp_gc, idx, bin)

      i = gcoff + gcinc
      @inbounds while i <= stop_i
         old_start = i - gcinc
         new_start = old_start + gcsize
         for offset in 0:(gcinc - 1)
            gc_cnt -= ifelse(_is_gc_2bit(raw_getindex(seq, old_start + offset)), 1, 0)
            gc_cnt += ifelse(_is_gc_2bit(raw_getindex(seq, new_start + offset)), 1, 0)
         end
         bin = Int(div(gc_cnt / inv_gcsize, gcwidth)) + 1
         gcfore[bin] += 1.0
         idx = _push_gc_temp!(temp_gc, idx, bin)
         i += gcinc
      end
   else
      @inbounds for i in gcoff:gcinc:stop_i
         gc_cnt = 0
         for pos in i:(i + gcsize - 1)
            gc_cnt += ifelse(_is_gc_2bit(raw_getindex(seq, pos)), 1, 0)
         end
         bin = Int(div(gc_cnt / inv_gcsize, gcwidth)) + 1
         gcfore[bin] += 1.0
         idx = _push_gc_temp!(temp_gc, idx, bin)
      end
   end

   return idx
end

function _gc_count!( mod::JointBiasMod, seq::BioSequence{A}, start_idx::Int ) where A <: BioSequences.Alphabet
   # GC fragment bias
   idx = start_idx
   gcsize = mod.gcsize
   gcwidth = mod.gcwidth
   gcfore = mod.gcfore
   temp_gc = mod.temp.gc
   inv_gcsize = Float64(gcsize)
   @inbounds for i in mod.gcoffset:mod.gcincrement:(length(seq)-gcsize)
      gc_cnt = 0
      for pos in i:(i + gcsize - 1)
         nt = seq[pos]
         gc_cnt += ifelse(nt == DNA_G || nt == DNA_C, 1, 0)
      end
      bin = Int(div( gc_cnt / inv_gcsize, gcwidth )+1)
      gcfore[bin] += 1.0
      idx = _push_gc_temp!(temp_gc, idx, bin)
   end

   return idx
end

function gc_count!( mod::JointBiasMod, seq::BioSequence{A} ) where A <: BioSequences.Alphabet
   end_idx = _gc_count!( mod, seq, 1 )
   if end_idx <= length(mod.temp.gc)
      @inbounds mod.temp.gc[end_idx] = zero(Int32)
   end
   return mod.temp
end

function count!( mod::JointBiasMod, seq::BioSequence{A} ) where A <: BioSequences.Alphabet
   # 5 prime bias
   primer_count!( mod, seq )
   # GC fragment bias
   end_idx = _gc_count!( mod, seq, 1 )
   if end_idx <= length(mod.temp.gc)
      @inbounds mod.temp.gc[end_idx] = zero(Int32)
   end

   return mod.temp
end

function count!( mod::JointBiasMod, fwd::BioSequence{A}, rev::BioSequence{A} ) where A <: BioSequences.Alphabet
   # 5 prime bias
   primer_count!( mod, fwd )
   # GC fragment bias (accumulate both mates into mod.temp.gc without overwriting fwd)
   mid_idx = _gc_count!( mod, fwd, 1 )
   end_idx = _gc_count!( mod, rev, mid_idx )
   if end_idx <= length(mod.temp.gc)
      @inbounds mod.temp.gc[end_idx] = zero(Int32)
   end

   return mod.temp
end

@fastmath value!( mod::JointBiasMod, hept::UInt16 ) = mod.fore[hept+1] > 0.0 ? mod.back[hept+1] / mod.fore[hept+1] : 1.0
@fastmath value!( mod::JointBiasMod, bin::Int ) = (mod.gcfore[bin] > 0.0 && mod.gcback[bin] > 0.0) ? mod.gcback[bin] / mod.gcfore[bin] : 1.0


mutable struct JointBiasCounter <: ReadCounter
   count::Float64
   map::Dict{UInt16,Int32}
   gc::Vector{Int32}
   isadjusted::Bool

   function JointBiasCounter( gc_param = GCBiasParams(50, 0.05) )
      gc = zeros(Int32, Int(ceil(1.0 / gc_param.width)+1))
      new(0.0, Dict{UInt16,Int32}(), gc, false)
   end
end

Base.zero(::Type{JointBiasCounter}) = JointBiasCounter()

function Base.push!( cnt::JointBiasCounter, temp::JointBiasTemp )
   for i in 1:length(temp.primer)
      push!( cnt, temp.primer[i] )
   end
   for i in 1:length(temp.gc)
      if temp.gc[i] > 0
         push!( cnt, temp.gc[i] )
      else
         return cnt
      end
   end
   cnt
end

function Base.push!( cnt::JointBiasCounter, bin::Int32 )
   cnt.gc[bin] += one(Int32)
   cnt
end

function Base.push!( cnt::JointBiasCounter, hept::UInt16 )
   if haskey( cnt.map, hept )
      cnt.map[hept] += 1
   else
      cnt.map[hept] = 1
   end
   cnt
end

adjust!( cnt::JointBiasCounter, mod::JointBiasMod ) = !cnt.isadjusted ? primer_adjust!( cnt, mod ) : gc_adjust!( cnt, mod )

@inline function primer_adjust!( cnt::JointBiasCounter, mod::JointBiasMod )
   for (k,v) in cnt.map
      cnt.count += value!( mod, k ) * v
   end
   cnt.count /= mod.forenpos
   cnt.isadjusted = true
end

@inline function gc_adjust!( cnt::JointBiasCounter, mod::JointBiasMod )
   weight = 0.0
   gcsum = Int32(0)
   gc = cnt.gc
   @inbounds for i in 1:length(gc)
      c = gc[i]
      if c > 0
         weight += value!( mod, i ) * c
         gcsum += c
      end
   end
   if gcsum > 0
      weight /= gcsum
      cnt.count *= weight
   end
   cnt
end
