
const SGKmer{K} = BioSequences.DNAMer{K}

Base.write(io::IO, k::BioSequences.DNAMer) = write(io, reinterpret(UInt64, k))
Base.read(io::IO, ::Type{T}) where {T <: BioSequences.DNAMer} = convert(T, Base.read(io, UInt64))

# Conversion
# ----------

kmer_index( kmer::BioSequences.DNAMer{K} ) where {K} = Int(reinterpret(UInt64, kmer)) + 1

kmer_index( seq::BioSequence ) = Int(kmer_index_trailing( seq )) + 1

# calculate kmer index directly
@inline function kmer_index_trailing( ::Type{T}, seq ) where T <: Unsigned
   x     = T(0)
   for nt in seq
      ntint = convert(UInt8, trailing_zeros(nt))
      if ntint > 0x03 || isambiguous(nt)
         return zero(T)
      else
         x = x << 2 | ntint
      end
   end
   x
end

kmer_index_trailing( seq ) = kmer_index_trailing(UInt64, seq )

# Zero-allocation kmer index over a subsequence seq[start:stop]
@inline function kmer_index_trailing_sub( ::Type{T}, seq::LongSequence{DNAAlphabet{2}}, start::Int, stop::Int ) where T <: Unsigned
   x = zero(T)
   @inbounds for i in start:stop
      x = (x << 2) | convert(T, raw_getindex(seq, i))
   end
   x
end

@inline function kmer_index_trailing_sub( ::Type{T}, seq::BioSequence, start::Int, stop::Int ) where T <: Unsigned
   x = zero(T)
   @inbounds for i in start:stop
      nt = seq[i]
      ntint = convert(UInt8, trailing_zeros(nt))
      if ntint > 0x03 || isambiguous(nt)
         return zero(T)
      else
         x = (x << 2) | convert(T, ntint)
      end
   end
   x
end

@inline kmer_index_sub( seq::BioSequence, start::Int, stop::Int ) = Int(kmer_index_trailing_sub(UInt64, seq, start, stop)) + 1

kmer(seq::SGSequence) = sgkmer(seq)
sgkmer(seq::String) = SGKmer(seq)
sgkmer(seq::SGSequence) = SGKmer(seq)
