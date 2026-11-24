using WaveletMatrices

for n in (2, 4)
    @eval begin
        bitsof(::Type{DNAAlphabet{$n}}) = $n
        bitsof(::Type{RNAAlphabet{$n}}) = $n
    end
end

@inline function raw_getindex(seq::BioSequence{A}, i::Integer) where A
    j = BioSequences.bitindex(seq, i)
    @inbounds return convert(UInt8, (seq.data[BioSequences.index(j)] >> BioSequences.offset(j)) & BioSequences.bitmask(bitsof(A)))
end

## This function is necessary as a result of changes made in Bio.jl as of Julia v0.5
## Since all nucleotides are converted to UInt8 as one hot 4-bit encodings by default
## we can override FMIndexes.sa_range function for those sequences specifically
## to re-encode as DNAAlphabet{2}
@inline function FMIndexes.sa_range(query::BioSequence{A}, index::FMIndex, init_range::UnitRange{Int}) where A <: BioSequences.Alphabet
    sp, ep = init_range.start, init_range.stop
    # backward search
    i = length(query)
    sentinel = index.sentinel
    count = index.count
    bwt = index.bwt
    @inbounds while sp ≤ ep && i ≥ 1
        char = raw_getindex(query, i)
        c = count[char+1]
        sp = c + WaveletMatrices.rank(char, bwt, (sp > sentinel ? sp - 1 : sp) - 1) + 1
        ep = c + WaveletMatrices.rank(char, bwt, (ep > sentinel ? ep - 1 : ep))
        i -= 1
    end
    return sp:ep
end

# Zero-allocation backward search over a subsequence query[qstart:qstop]
@inline function sa_range_sub(query::LongSequence{DNAAlphabet{2}}, qstart::Int, qstop::Int, index::FMIndex)
    sp, ep = 1, length(index) + 1
    i = qstop
    sentinel = index.sentinel
    count = index.count
    bwt = index.bwt
    @inbounds while sp ≤ ep && i ≥ qstart
        char = raw_getindex(query, i)
        c = count[char+1]
        sp = c + WaveletMatrices.rank(char, bwt, (sp > sentinel ? sp - 1 : sp) - 1) + 1
        ep = c + WaveletMatrices.rank(char, bwt, (ep > sentinel ? ep - 1 : ep))
        i -= 1
    end
    return sp:ep
end

# Zero-allocation backward search over reverse_complement(query[qstart:qstop])
@inline function sa_range_sub_rc(query::LongSequence{DNAAlphabet{2}}, qstart::Int, qstop::Int, index::FMIndex)
    sp, ep = 1, length(index) + 1
    i = qstart
    sentinel = index.sentinel
    count = index.count
    bwt = index.bwt
    @inbounds while sp ≤ ep && i ≤ qstop
        char = raw_getindex(query, i) ⊻ 0x03
        c = count[char+1]
        sp = c + WaveletMatrices.rank(char, bwt, (sp > sentinel ? sp - 1 : sp) - 1) + 1
        ep = c + WaveletMatrices.rank(char, bwt, (ep > sentinel ? ep - 1 : ep))
        i += 1
    end
    return sp:ep
end
