
const _DNA2_SANITIZE = let tbl = fill(UInt8('A'), 256)
   for (c, v) in (('A','A'), ('a','A'), ('C','C'), ('c','C'),
                  ('G','G'), ('g','G'), ('T','T'), ('t','T'))
      tbl[Int(UInt8(c)) + 1] = UInt8(v)
   end
   tbl
end

mutable struct FASTQRecord
   sequence::LongSequence{DNAAlphabet{2}}
   quality::Vector{UInt8}
   raw::FASTQ.Record

   FASTQRecord() = new( LongSequence{DNAAlphabet{2}}(),
                        Vector{UInt8}(), FASTQ.Record() )
end

Base.read!( parser::FASTQ.Reader, rec::FASTQRecord ) = read!( parser, rec.raw )

function Base.fill!( rec::FASTQRecord, offset=33 )
   raw = rec.raw
   data = raw.data
   srange = raw.sequence
   slen = length(srange)
   soff = first(srange) - 1
   sanitize = _DNA2_SANITIZE
   @inbounds for i in 1:slen
      data[soff + i] = sanitize[Int(data[soff + i]) + 1]
   end
   seq = rec.sequence
   resize!(seq, slen)
   copyto!(seq, 1, data, first(srange), slen)

   qrange = raw.quality
   qlen = length(qrange)
   qual = rec.quality
   resize!(qual, qlen)
   qoff = first(qrange) - 1
   off_u8 = UInt8(offset)
   @inbounds for i in 1:qlen
      qual[i] = data[qoff + i] - off_u8
   end
   rec
end

Base.string( rec::FASTQRecord ) = String(rec.raw.data[rec.raw.identifier])

function reverse_complement!(rec::FASTQRecord)
   BioSequences.reverse_complement!(rec.sequence)
   reverse!(rec.quality)
end
