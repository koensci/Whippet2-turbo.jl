using Test
using Serialization
using Pkg

import BioSymbols
import BioSequences
import BioAlignments
import XAM
import FASTX
using BioSymbols
using BioSequences
using BioAlignments
using GenomicFeatures
using XAM
using FASTX
using BufferedStreams
using DataStructures
using Distributions
using FMIndexes
using IntArrays
using InteractiveUtils
using IntervalTrees
using Libz
using LinearAlgebra
using Nullables
using Printf
using Random
using CSV
using DataFrames

include("../src/types.jl")
include("../src/timer.jl")
include("../src/sgkmer.jl")
include("../src/fmindex_patch.jl")
include("../src/bam.jl")
include("../src/refset.jl")
include("../src/graph.jl")
include("../src/bias.jl")
include("../src/edges.jl")
include("../src/index.jl")
include("../src/record.jl")
include("../src/align.jl")
include("../src/quant.jl")
include("../src/reads.jl")
include("../src/paired.jl")
include("../src/events.jl")
include("../src/io.jl")
include("../src/diff.jl")
include("../src/diff_glm.jl")

@testset "SG Kmers & Seq" begin
   @test sgkmer(dna"ATG") == BioSequences.DNAMer(dna"ATG")
   @test sgkmer("ATG") == BioSequences.DNAMer(dna"ATG")
   @test isa(sgkmer(dna"ATG"), SGKmer{3})
   @test kmer_index( sgkmer(dna"ATG") ) == kmer_index( BioSequences.DNAMer(dna"ATG") )
   @test Int(kmer_index_trailing( dna"ATG" ))+1 == kmer_index( BioSequences.DNAMer(dna"ATG") )
   @test Int(kmer_index_trailing( dna"ATGR" )) == 0
   seq = dna""
   for i in 1:32
      seq = i % 2 == 0 ? seq * dna"C" : seq * dna"T"
      curkmer = sgkmer(seq)
      @test isa( curkmer, SGKmer{i} )
   end
end

@testset "Bias Models" begin
   # default 'empty' model
   def = zero(DefaultCounter)
   @test get(def) == 0.0
   def = one(DefaultCounter)
   @test get(def) == 1.0
   push!( def, 1.0 )
   @test get(def) == 2.0
   gc_adjust!( def, DefaultBiasMod() )
   primer_adjust!( def, DefaultBiasMod() )

   function randdna(n)
      return BioSequence{DNAAlphabet{2}}(rand([DNA_A, DNA_C, DNA_G, DNA_T], n))
   end

   function test_uniform_bias!( mod::B ) where B <: BiasModel
      for i in 1:5000000
         primer_count!( mod, randdna(36) )
      end
      primer_normalize!( mod )
      for i in 1:length(mod.fore)
         @test 0.85 < mod.back[i] / mod.fore[i] < 1.15
      end
   end

   # primer count model
   mod = JointBiasMod( 5 )
   test_uniform_bias!( mod )
   cnt = JointBiasCounter()
   @test cnt.count == 0.0
   hept = kmer_index_trailing(UInt16, dna"ATGAC")
   push!( cnt, hept )
   push!( cnt, hept )
   adjust!( cnt, mod )
   @test get(cnt) == value!( mod, hept )

   # gc content bias
   seq = dna"AAAAAGCGCG" ^ 5
   egc = ExpectedGC( seq )
   @test length(egc) == 20
   @test egc[ Int(div(gc_content(seq), 0.05)+1) ] == 1.0

   # Regression tests for GC bias lookup (gcfore/gcback) and paired-end GC bin accumulation
   gcmod = JointBiasMod()
   # Distinct GC contents: low-GC fwd (~20% GC -> bin 5) and high-GC rev (~80% GC -> bin 17)
   fwd_seq = LongSequence{DNAAlphabet{2}}("AAAAG" ^ 20)
   rev_seq = LongSequence{DNAAlphabet{2}}("GGGGCT" ^ 16 * "GGGG")
   tmp = count!( gcmod, fwd_seq, rev_seq )
   gccnt = JointBiasCounter()
   push!( gccnt, tmp )
   fwd_bin = Int(div(0.2, 0.05)) + 1
   rev_bin = Int(div(gc_content(rev_seq[5:54]), 0.05)) + 1
   @test gccnt.gc[fwd_bin] > 0
   @test gccnt.gc[rev_bin] > 0
   # Ensure value!(gcmod, bin::Int) reads gcback/gcfore rather than primer fore/back
   gcmod.fore .= 1.0
   gcmod.back .= 100.0
   gcmod.gcfore .= 0.0
   gcmod.gcback .= 0.0
   gcmod.gcfore[fwd_bin] = 0.25
   gcmod.gcback[fwd_bin] = 0.5
   @test value!( gcmod, fwd_bin ) ≈ 2.0
   @test value!( gcmod, rev_bin ) == 1.0
   @test add_to_back!( gcmod.gcback, egc, 1.0 ) === gcmod.gcback

   # Regression test for SRA FASTQ records with +SRR... line-3 text containing N/n bases
   sra_fq = IOBuffer("@SRR8615273.1 1 length=10\nNCGTnCGTAN\n+SRR8615273.1 1 length=10\n#III#IIII#\n")
   sra_reader = FASTQ.Reader( BufferedInputStream(sra_fq), fill_ambiguous=DNA_A )
   sra_rec = FASTQRecord()
   read!( sra_reader, sra_rec )
   fill!( sra_rec, 33 )
   @test sra_rec.sequence == LongSequence{DNAAlphabet{2}}("ACGTACGTAA")
   @test sra_rec.quality == UInt8[2, 40, 40, 40, 2, 40, 40, 40, 40, 2]
end

@testset "Splice Graphs" begin
   gtf = IOBuffer("# gtf file test
chr0\tTEST\texon\t6\t20\t.\t+\t.\tgene_id \"one\"; transcript_id \"def\";
chr0\tTEST\texon\t31\t40\t.\t+\t.\tgene_id \"one\"; transcript_id \"def\";
chr0\tTEST\texon\t54\t62\t.\t+\t.\tgene_id \"one\"; transcript_id \"def\";
chr0\tTEST\texon\t76\t85\t.\t+\t.\tgene_id \"one\"; transcript_id \"def\";
chr0\tTEST\texon\t6\t40\t.\t+\t.\tgene_id \"one\"; transcript_id \"int1_alt3\";
chr0\tTEST\texon\t51\t62\t.\t+\t.\tgene_id \"one\"; transcript_id \"int1_alt3\";
chr0\tTEST\texon\t76\t85\t.\t+\t.\tgene_id \"one\"; transcript_id \"int1_alt3\";
chr0\tTEST\texon\t6\t20\t.\t+\t.\tgene_id \"one\"; transcript_id \"apa_alt5\";
chr0\tTEST\texon\t31\t40\t.\t+\t.\tgene_id \"one\"; transcript_id \"apa_alt5\";
chr0\tTEST\texon\t54\t65\t.\t+\t.\tgene_id \"one\"; transcript_id \"apa_alt5\";
chr0\tTEST\texon\t76\t90\t.\t+\t.\tgene_id \"one\"; transcript_id \"apa_alt5\";
chr0\tTEST\texon\t11\t20\t.\t+\t.\tgene_id \"single\"; transcript_id \"ex1_single\";
chr0\tTEST\texon\t11\t20\t.\t-\t.\tgene_id \"single_rev\"; transcript_id \"single_rev\";
chr0\tTEST\texon\t11\t20\t.\t-\t.\tgene_id \"kissing\"; transcript_id \"def_kiss\";
chr0\tTEST\texon\t21\t30\t.\t-\t.\tgene_id \"kissing\"; transcript_id \"def_kiss\";
chr0\tTEST\texon\t11\t30\t.\t-\t.\tgene_id \"kissing\"; transcript_id \"ret_kiss\";
chr0\tTEST\texon\t11\t30\t.\t-\t.\tgene_id \"spl_neg\"; transcript_id \"spl_negative\";
chr0\tTEST\texon\t92\t97\t.\t-\t.\tgene_id \"spl_neg\"; transcript_id \"spl_negative\";
")


#1         11        21        31        41        51        61        71        81        91         101
#|.........|.........|.........|.........|.........|.........|.........|.........|.........|..........|
#     [             ]----------[        ]-------------[       ]-------------[        ]       + def
#     [             .          .        ]----------[  .       ]-------------[        ]       + int1_alt3
#     [             ]----------[        ]-------------[       .  ]----------[        .    ]  + apa_alt5
#
#          [        ] + ex1_single
#
#          [        ] - single_rev
#
#          [        ][        ] - def_kiss
#          [                  ] - ret_kiss
#          [                  ]-------------------------------------------------------------[     ] - spl_negative
#            
#

   gtfref  = load_gtf( gtf )

   @testset "Gene Annotations" begin

      @test gtfref["one"].don == (20,40,62,65)
      @test gtfref["one"].acc == (31,51,54,76)
      @test gtfref["one"].noveldon == tuple()
      @test gtfref["one"].novelacc == tuple()

   end

   gtf.ptr = 1
   bamreadr = open(BAM.Reader, "test.sort.bam", index="test.sort.bam.bai" )
   supgtfref = load_gtf( gtf, usebam=true, bamreader=Nullable(bamreadr), bamreads=2 )

   @testset "BAM Supplemented Annotation" begin

      @test supgtfref["one"].don == (20,40,45,62,65)
      @test supgtfref["one"].acc == (31,51,54,76)
      @test supgtfref["one"].noveldon == tuple(45)
      @test supgtfref["one"].novelacc == tuple()

   end

                               # fwd     rev
   buffer1   = dna"AAAAA"      # 1-5     96-100
   utr5      = dna"TTATT"      # 6-10    91-95
   exon1     = dna"GCGGATTACA" # 11-20   81-90
   int1      = dna"TTTTTTTTTT" # 21-30   71-80
   exon2     = dna"GCATTAGAAG" # 31-40   61-70
   int2      = dna"GGGGGGGGGG" # 41-50   51-60
   exon3alt3 = dna"CCT"        # 51-53   48-50
   exon3def  = dna"CTATGCTAG"  # 54-62   39-47
   exon3alt5 = dna"TTC"        # 63-65   36-38
   int3      = dna"CCCCCCCCCC" # 66-75   26-35
   exon4     = dna"TTAGACAAGA" # 76-85   16-25
   apa       = dna"AATAA"      # 86-90   11-15
   buffer2   = dna"AAAAAAAAAA" # 91-100  1-10

   fwd = buffer1 * utr5 * exon1 * int1 *
         exon2 * int2 *
         exon3alt3 * exon3def * exon3alt5 * int3 *
         exon4 * apa * buffer2
   rev = reverse_complement(fwd)

   genome = fwd * rev

   expected_one = utr5 * exon1 *
               int1 *
               exon2 *
               exon3alt3 *
               exon3def *
               exon3alt5 *
               exon4 *
               apa

   expected_sin = dna"GCGGATTACA"
   expected_kis = dna"AAAAAAAAAATGTAATCCGC"

   kmer_size = 2 # good test size

   graph_one = SpliceGraph( gtfref["one"], genome, kmer_size )
   graph_sin = SpliceGraph( gtfref["single"], genome, kmer_size )
   graph_rev = SpliceGraph( gtfref["single_rev"], genome, kmer_size )
   graph_kis = SpliceGraph( gtfref["kissing"], genome, kmer_size )

   @testset "Graph Building" begin
      @test graph_one.seq == expected_one
      @test graph_sin.seq == expected_sin
      @test graph_kis.seq == expected_kis

      @test length(graph_one.annopath) == length(gtfref["one"].reftx)
      @test graph_one.annopath[1] == BitSet([1,3,5,7]) # path of def
      @test graph_one.annopath[2] == BitSet([1,2,3,4,5,7]) # path of int1_alt3
      @test graph_one.annopath[3] == BitSet([1,3,5,6,7,8]) # path of apa_alt5
   end

   # Build Index (from index.jl)
   xcript  = dna""
   xoffset = Vector{UInt64}()
   xgenes  = Vector{GeneName}()
   xinfo   = Vector{GeneInfo}()
   xlength = Vector{Float64}()
   xgraph  = Vector{SpliceGraph}()

   runoffset = 0

   for g in keys(gtfref)
      println(stderr, g)
      curgraph = SpliceGraph( gtfref[g], genome, kmer_size )
      xcript  *= curgraph.seq
      push!(xgraph, curgraph)
      push!(xgenes, g)
      push!(xinfo, gtfref[g].info )
      push!(xlength, gtfref[g].length )
      push!(xoffset, runoffset)
      runoffset += length(curgraph.seq)
   end

   fm = FMIndex(twobit_enc(xcript), 4, r=1, program=:SuffixArrays, mmap=true)

   edges = build_edges( xgraph, kmer_size )

   #println(edges.left)
   #println(edges.right)
   @time ic = interval_index( xgraph, xinfo, min_intron_size=1 )

   lib = GraphLib( xoffset, xgenes, xinfo, xlength, xgraph, ic, edges, fm, true, kmer_size, true, v"1.7" )
   @test checkversion(lib, v"1.7", v"1.7")

   @testset "Kmer Edges" begin
      left  = [dna"CA", dna"AG", dna"AG", dna"TC", dna"AA", dna"TT"]
      right = [dna"GC", dna"CC", dna"CT", dna"TT", dna"TG", dna"AA"]
      lkmer = map( x->kmer_index(SGKmer(x)), left )
      rkmer = map( x->kmer_index(SGKmer(x)), right )
      #println(lkmer[1])
      #println(rkmer[1])
      for i in 1:4^kmer_size
         if i in lkmer
            @test isassigned(edges.left, i)
            @test typeof(edges.left[i]) == Vector{SGNode}
            @test issorted(edges.left[i], lt=sortlt)
         else
            @test !isassigned(edges.left, i)
         end
         if i in rkmer
            @test isassigned(edges.right, i)
            @test typeof(edges.right[i]) == Vector{SGNode}
            @test issorted(edges.right[i], lt=sortlt)
         else
            @test !isassigned(edges.right, i)
         end
      end
      exon1_lind = lkmer[1]
      exon2_rind = rkmer[1]
      @test intersect( edges.left[exon1_lind], edges.right[exon2_rind] ) == edges.right[exon2_rind]
      @test edges.left[exon1_lind] ∩ edges.right[exon2_rind] == edges.right[exon2_rind]
   end

   @testset "Saving and Loading Index" begin
      println(stderr, "Saving test index...")
      println(stderr, lib)
      @timer open("test_index.jls", "w+") do io
         serialize( io, lib )
      end
      println(stderr, "Loading test index...")
      @timer lib = open(deserialize, "test_index.jls")
      println(stderr, lib)
   end

   @testset "Building IntervalCollection for BAM" begin
      ic = interval_index(lib.graphs, lib.info, min_intron_size=1)
      @test length(collect(ic)) == 17 # total nodes of len >= 1
      @test sum(map(x->isexonic(x.metadata), collect(ic))) == 14 # two are introns
   end

   # store general data structures outside of testset
   score_range = 0.05
   param = AlignParam( 1, 2, 4, 4, 4, 5, 1, 2, 1000, score_range, 0.7,
                          false, false, true, false, true )
   println(map( x->length(x.annoname), lib.graphs ))
   gquant = GraphLibQuant{SGAlignSingle,DefaultCounter}( lib )

   @testset "Alignment, SAM Format, Equivalence Classes" begin
      # reads
      fastq_str = "@1S10M%11,20%exon1
NGCGGATTACA
+
#BBBBBBBBBB
@1S9M%54,62%exon3def
NCTATGCTAG
+
#BBBBBBBBB
@1S15M%51,65%alt3-exon3-alt5
NCCTCTATGCTAGTTC
+
#BBBBBBBBBBBBBBB
@11M10N10M%10,40%exon1-exon2
NGCGGATTACAGCATTAGAAG
+
#BBBBBBBBBBBBBBBBBBBB
@5M10N6M%16,36%exon1trunc-exon2trunc
TTACAGCATTN
+
BBBBBBBBBB#
@10M33N9M1S%11,62%exon1-exon3def
GCGGATTACACTATGCTAGN
+
BBBBBBBBBBBBBBBBBBB#
@10M33N9M%11,62%exon1-exon3def:rc
CTAGCATAGTGTAATCCGC
+
BBBBBBBBBBBBBBBBBBB
@10M55N10M1S%11,85%exon1-exon4full
GCGGATTACATTAGACAAGAN
+
BBBBBBBBBBBBBBBBBBBB#
@11M55N4M%10,78%exon1-exon4_4bp
NGCGGATTACATTAG
+
#IIIIIIIIIIIIII
@2M55N10M%19,85%exon1_2bp-exon4:rc
TCTTGTCTAATG
+
IIIIIIIIIIII
"
      fastq = IOBuffer(fastq_str)
      DNASeqType = BioSequences.BioSequence{BioSequences.DNAAlphabet{2}}
      fqparse = FASTQ.Reader( BufferedInputStream(fastq), fill_ambiguous=DNA_C )
      reads  = allocate_fastq_records( 10 )
      read_chunk!( reads, fqparse )
      sambuf = BufferedOutputStream( open("test_out.sam", "w") )
      write_sam_header( sambuf, lib )

      aligns = Vector{SGAlignment}()

      @test length(reads) == 10
      for r in reads
         fill!( r, 33 )
         println(stderr, r)
         readname = string(r)

         align = ungapped_align( param, lib, r )
         #println(stderr, align)

         flush(stderr)
         @test !isnull( align )
         @test length( align.value ) >= 1
         @test all(map( x->x.isvalid, align.value))
         scores = map( x->identity(x, length(r.sequence)), align.value )
         @test maximum(scores) - minimum(scores) <= score_range

         if length(something(findfirst(":rc", readname), 0:-1)) > 0
            @test align.value[1].strand == false
         else
            @test align.value[1].strand == true
         end

         best_ind = argmax(scores)

         ex_num = length(split(readname, '-', keepempty=false))
         @test length(align.value[best_ind].path) == ex_num

         println(stderr, align.value[best_ind])
         count!( gquant, align.value[best_ind], 1.0 )

         curgene   = align.value[best_ind].path[1].gene
         firstnode = align.value[best_ind].path[1].node
         lastnode  = align.value[best_ind].path[end].node
         curgraph  = lib.graphs[ curgene ]

         cigar,positions = split(readname, '%', keepempty=false)[1:2]
         first,last   = split(positions, ',', keepempty=false) |> y->map(x->parse(Int,x), y)
         # Test SAM Format
         @test cigar == cigar_string( align.value[best_ind], curgraph, true, length(r.sequence) )[1]
         test_cigar,endpos = cigar_string( align.value[best_ind], curgraph, true, length(r.sequence) )
         #println(stderr, "cigar = $test_cigar")
         @test first == Int(curgraph.nodecoord[firstnode] + (align.value[best_ind].offset - curgraph.nodeoffset[firstnode]))
         #@test last  == Int(curgraph.nodecoord[lastnode] + (endpos - curgraph.nodeoffset[lastnode]))
         # test readlength = number of M and S entries in cigar
         # test SAM offset is correct for both '+' and '-' genes.
         # test that cigar reversal works for '-' strand genes.
         push!(aligns, align.value[best_ind])

         write_sam( sambuf, r, align.value[best_ind], lib, qualoffset=33 )
      end
      close(sambuf)

      align_quant = deepcopy(gquant)
      bam_quant = GraphLibQuant{SGAlignSingle,DefaultCounter}( lib )

      @testset "Aligning SAM/BAM files" begin
         parser = open(SAM.Reader, "test_out.sam")
         data = AlignData() 
         overlaps  = Dict{NodeInt, Int}()
         mod = DefaultBiasMod()
         for r in parser
            path, valid = align_bam( lib.coords, lib.info, data, overlaps, r, true, 1 )
            println(stderr, path)
            cur_sam = popfirst!( aligns )
            @test valid
            @test path == cur_sam.path
            @test hash(path) == hash(cur_sam.path)
            println(stderr, path)
            count!( bam_quant, path, 0.0 )
            seq_fwd = isstrandfwd(r) ? sequence(r) : reverse_complement(sequence(r))
            biasval = count!( mod, seq_fwd )
            count!( bam_quant, path, biasval )
         end
         @test length(aligns) == 0
         close(parser)
         @test ispaired_bamfile( "test.sort.bam" ) == false
      end

      @testset "Aberrant Reads" begin

         for i in 1:99999
            enc = encode_aberrant(i, i)
            root, node, pos = decode_aberrant( enc )
            @test i == root
            @test i == pos
         end
         for i in 1:0.5:99999
            enc = encode_aberrant(i, 12345)
            root, node, pos = decode_aberrant( enc )
            @test i == node
            @test pos == 12345
         end

      end

      @testset "SGAlignContainer Hashing" begin

         @test SGAlignNode(1,1,zero(SGAlignScore)) == SGAlignNode(1,1,one(SGAlignScore))
         @test hash( SGAlignNode(1,1,zero(SGAlignScore)) ) == hash( SGAlignNode(1,1,one(SGAlignScore)) )

         @test SGAlignNode[SGAlignNode(1,1,zero(SGAlignScore))] == SGAlignNode[SGAlignNode(1,1,one(SGAlignScore))]
         @test hash( SGAlignNode[SGAlignNode(1,1,zero(SGAlignScore))] ) == hash( SGAlignNode[SGAlignNode(1,1,one(SGAlignScore))] )

         @test SGAlignSingle(SGAlignNode[SGAlignNode(1,1,zero(SGAlignScore))]) == SGAlignSingle(SGAlignNode[SGAlignNode(1,1,one(SGAlignScore))])
         @test hash( SGAlignSingle(SGAlignNode[SGAlignNode(1,1,zero(SGAlignScore))]) ) == hash( SGAlignSingle(SGAlignNode[SGAlignNode(1,1,one(SGAlignScore))]) )
     
         path = SGNodeMeta{NodeInt,Nothing}[SGNodeMeta(NodeInt(1),NodeInt(1),nothing)]

         @test SGAlignPaired( path ) == SGAlignPaired( path, path )
         @test hash(SGAlignSingle(path)) == hash(SGAlignPaired(path))

         a = map( x->SGAlignSingle(SGAlignNode[SGAlignNode(x,x,zero(SGAlignScore))]), collect(1:5))
         ab = map( x->SGAlignSingle(SGAlignNode[SGAlignNode(x,x,one(SGAlignScore))]), collect(1:5))
         b = map( x->SGAlignPaired(SGAlignNode[SGAlignNode(x,x,zero(SGAlignScore))], SGAlignNode[SGAlignNode(x,x,zero(SGAlignScore))]), collect(1:5))
         bb = map( x->SGAlignPaired(SGAlignNode[SGAlignNode(x,x,one(SGAlignScore))], SGAlignNode[SGAlignNode(x,x,one(SGAlignScore))]), collect(1:5))

         da = Dict{SGAlignSingle,Int}()
         db = Dict{SGAlignPaired,Int}()

         function settwice!( d, a, ab )
            for i in a
               d[i] = 1
            end
            for i in ab
               d[i] = 2
            end
            for k in keys(d)
               @test d[k] == 2
            end
         end
         settwice!( da, a, ab )
         settwice!( db, b, bb )
         @test length(keys(da)) == length(a)
         @test length(keys(db)) == length(b)

      end

      @testset "SGAlignPath Overlap" begin

         single     = SGAlignNode[SGAlignNode(1,1,one(SGAlignScore))]
         single_two = SGAlignNode[SGAlignNode(1,2,one(SGAlignScore))]
         two_three  = SGAlignNode[SGAlignNode(1,2,one(SGAlignScore)), SGAlignNode(1,3,one(SGAlignScore))]
         larger     = SGAlignNode[SGAlignNode(1,1,one(SGAlignScore)), SGAlignNode(1,2,one(SGAlignScore)), SGAlignNode(1,3,one(SGAlignScore))]

         println(stderr, @which single in single)

         @test single in single
         @test !(single in single_two)
         @test single_two in two_three
         @test !(single in two_three)
         @test single in larger
         @test single_two in larger
         @test two_three in larger

      end

      @testset "Isoform Equivalence Classes" begin

         is = BitSet([1, 2, 5, 6])
         @test in( 1, 2, is ) == true
         @test in( 1, 3, is ) == false
         @test in( 2, 3, is ) == false
         @test in( 2, 4, is ) == false
         @test in( 2, 5, is ) == true
         @test in( 5, 6, is ) == true
         @test in( 1, 6, is ) == false
         for i in [1,2,5,6]
            @test (i in is) == true
         end
         path = SGAlignNode[SGAlignNode(0x00000002, 0x00000001, SGAlignScore(0x02, 0x00, 0.0)), SGAlignNode(0x00000002, 0x00000007, SGAlignScore(0x0a, 0x00, 0.0))]
         @test (path in is) == false
         path = SGAlignNode[SGAlignNode(0x00000002, 0x00000001, SGAlignScore(0x02, 0x00, 0.0)), SGAlignNode(0x00000002, 0x00000002, SGAlignScore(0x0a, 0x00, 0.0))]
         @test (path in is) == true
         path = SGAlignNode[SGAlignNode(0x00000002, 0x00000001, SGAlignScore(0x02, 0x00, 0.0)), SGAlignNode(0x00000002, 0x00000002, SGAlignScore(0x0a, 0x00, 0.0)), SGAlignNode(0x00000002, 0x00000003, SGAlignScore(0x0a, 0x00, 0.0))]
         @test (path in is) == false
         path = SGAlignNode[SGAlignNode(0x00000002, 0x00000001, SGAlignScore(0x02, 0x00, 0.0)), SGAlignNode(0x00000002, 0x00000002, SGAlignScore(0x0a, 0x00, 0.0)), SGAlignNode(0x00000002, 0x00000005, SGAlignScore(0x0a, 0x00, 0.0))]
         @test (path in is) == true

         @test length(gquant.tpm) == 8
         @test length(gquant.count) == 8
         @test length(gquant.length) == 8
         @test gquant.geneidx == [0, 2, 5, 6, 7] 
         build_equivalence_classes!( gquant, lib, assign_long=true )
         println(stderr, gquant)
         println(stderr, gquant.classes)

      end

      @testset "MultiMapping Equivalence Classes" begin

         multi  = MultiMapping{SGAlignSingle,DefaultCounter}()
         aligns = SGAlignment[SGAlignment(0x0000000e, 0x01, SGAlignNode[SGAlignNode(0x00000002, 0x00000001, SGAlignScore(0x02, 0x00, 0.0)), SGAlignNode(0x00000002, 0x00000007, SGAlignScore(0x0a, 0x00, 0.0))], false, true),
                              SGAlignment(0x00000001, 0x01, SGAlignNode[SGAlignNode(0x00000004, 0x00000001, SGAlignScore(0x02, 0x00, 0.0))], false, true)]

         push!( multi, aligns, 1.0, gquant, lib )
         println(stderr, multi.map)
         compats = collect(values(multi.map))
         @test compats[1].isdone == true
         @test sum(multi.iset) == 0

         multi  = MultiMapping{SGAlignSingle,DefaultCounter}()

         aligns = SGAlignment[SGAlignment(0x0000000e, 0x01, SGAlignNode[SGAlignNode(0x00000002, 0x00000001, SGAlignScore(0x02, 0x00, 0.0)), SGAlignNode(0x00000002, 0x00000002, SGAlignScore(0x0a, 0x00, 0.0))], false, true),
                              SGAlignment(0x00000001, 0x01, SGAlignNode[SGAlignNode(0x00000002, 0x00000001, SGAlignScore(0x02, 0x00, 0.0))], false, true)]

         push!( multi, aligns, 1.0, gquant, lib )
         println(stderr, multi.map)
         compats = collect(values(multi.map))
         @test compats[1].isdone == false
         @test [3,4,5] == compats[1].class
         @test sum(multi.iset) == 0

         # Assignment
         compats[1].count = 10.0
         prev_node = get(gquant.quant[2].node[1])
         println(stderr, "prev_quant = $(gquant.quant[2])")
         interv = IntervalTrees.Interval{NodeInt}( 1, 2 )
         prev_edge = get(get( gquant.quant[2].edge, interv, IntervalValue(0,0,zero(DefaultCounter)) ).value)

         println(gquant.quant[2].edge)

         assign_ambig!( gquant, lib, multi )
         println(gquant.quant[2].edge)

         println(stderr, "cur_quant = $(gquant.quant[2])")
         @test get(gquant.quant[2].node[1]) == prev_node + 7.5
         @test get(get(gquant.quant[2].edge, interv, IntervalValue(0,0,default(DefaultCounter)) ).value) == prev_edge + 2.5
      end

      @testset "TPM" begin
         multi  = MultiMapping{SGAlignSingle,DefaultCounter}()
         calculate_tpm!( gquant, readlen=20 )
         set_gene_tpm!( gquant, lib )
         orig = deepcopy(gquant.tpm)
         iter = gene_em!( gquant, multi, sig=1, readlen=20, maxit=10000 )
         @test orig != gquant.tpm
         @test iter < 5000
         @test sum(gquant.tpm) == 1000000
      end

      @testset "Output Files" begin
         output_junctions( "test.jnc.gz", lib, gquant )
         output_stats( "test.map.gz", lib, gquant, param, "test_index.jls", 10, 10, 0, 20, "vTEST" )
         output_exons( "test_index.exons.gz", lib )
      end

      function parse_edge( str::S ) where S <: AbstractString
         s = split( str, ['-',':'] )
         (parse(Int, String(s[1])), parse(Int, String(s[2])), parse(Float64, String(s[3])))
      end

      function psigraph_from_string( edgestr::String )
         edgespl  = split( edgestr, ',' )
         edges    = IntervalMap{Int,Float64}()
         psigraph = PsiGraph{BitSet}( Vector{Float64}(), Vector{Float64}(),
                                      Vector{Float64}(), Vector{BitSet}(), 1, 10 )
         for s in split( edgestr, ',' )
            f,l,v = parse_edge( s )
            edges[(f,l)] = v
            push!( psigraph.psi, 0.0 )
            push!( psigraph.length, 1.0 )
            push!( psigraph.count, v )
            push!( psigraph.nodes, BitSet([f,l]) )
         end
         psigraph
      end

      @testset "Graph Enumeration" begin
         #= Real Test based on node 23 (here annotated as 1) in ENSMUSG00000038685
         inclusion edges:
         1-2, 2-3, 3-4, 4-5, 5-6, 6-7, 7-8, 8-9, 9-10
         exclusion edges:
         1-3, 2-4, 3-6, 3-10, 4-6, 8-10
         expected paths:
         =#
         expected_paths = [BitSet([1,2,3,4,5,6,7,8,9,10]) # (full)
         BitSet([1,3,4,5,6,7,8,9,10]) # (1-3)
         BitSet([1,2,4,5,6,7,8,9,10]) # (2-4)
         BitSet([1,2,3,6,7,8,9,10])   # (3-6)
         BitSet([1,3,6,7,8,9,10])     # (3-6)
         BitSet([1,2,3,10])           # (3-10)
         BitSet([1,3,10])             # (3-10)
         BitSet([1,2,3,4,6,7,8,9,10]) # (4-6)
         BitSet([1,3,4,6,7,8,9,10])   # (4-6)
         BitSet([1,2,4,6,7,8,9,10])   # (4-6)
         BitSet([1,2,3,4,5,6,7,8,10]) # (8-10)
         BitSet([1,3,4,5,6,7,8,10])   # (8-10)
         BitSet([1,2,4,5,6,7,8,10])   # (8-10)
         BitSet([1,2,3,6,7,8,10])     # (8-10)
         BitSet([1,3,6,7,8,10])       # (8-10)
         BitSet([1,2,3,4,6,7,8,10])   # (8-10)
         BitSet([1,3,4,6,7,8,10])     # (8-10)
         BitSet([1,2,4,6,7,8,10])]    # (8-10)
         edgestr = "1-2:616.0,1-3:4.0,2-3:569.0,2-4:20.0,3-4:629.0,3-6:1.0,3-10:3.0,4-5:1.0,4-6:664.0,5-6:5.0,6-7:789.0,7-8:790.0,8-9:2.0,8-10:606.0,9-10:15.0"
         psigraph = psigraph_from_string( edgestr )
         res = reduce_graph( psigraph )
         for p in expected_paths
            @test p in res.nodes
         end
         @test length(expected_paths) == length(res.nodes)

         expected_paths = [BitSet([6,14,15,16])
         BitSet([6,14,17])
         BitSet([6,17])
         BitSet([7,8,9,10,11,12,13,14,15,16])
         BitSet([7,8,9,10,11,12,13,14,17])
         BitSet([7,8,9,10,11,12,13,17])]
         nrf2test = "6-14:129.0,6-17:3.0,7-8:2.0,8-9:5.0,9-10:8.0,10-11:4.0,11-12:26.0,12-13:22.0,13-14:36.0,13-17:11.0,14-15:3.0,14-17:382.0,15-16:7.0"
         psigraph = psigraph_from_string( nrf2test )

         res = reduce_graph( psigraph )
         for p in expected_paths
            @test p in res.nodes
         end
         @test length(expected_paths) == length(res.nodes)

      end

      @testset "Event Building" begin
         out = IOBuffer(read=true,write=true)
         bs  = BufferedOutputStream(out)
         output_psi_header( bs )
         for g in 1:length(lib.graphs)
            name = lib.names[g]
            chr  = lib.info[g].name
            strand = lib.info[g].strand ? '+' : '-'
            _process_events( bs, lib.graphs[g], gquant.quant[g], (name,chr,strand), readlen=20 )
         end
         flush(bs)
         seek(out,0)
         for l in eachline(out)
            println(stderr,l)
            spl = split( l, '\t' )
            @test length(spl) == 14
         end
         println(stderr,"")

         psi_vals = []
         conf_vals = []
         out = IOBuffer(read=true,write=true)
         bs  = BufferedOutputStream(out)
         output_psi_header( bs )
         for g in 1:length(lib.graphs)
            name = lib.names[g]
            chr  = lib.info[g].name
            strand = lib.info[g].strand ? '+' : '-'
            _process_events( bs, lib.graphs[g], align_quant.quant[g], (name,chr,strand), readlen=20 )
         end
         flush(bs)
         seek(out,0)
         for l in eachline(out)
            println(stderr,l)
            spl = split( l, '\t' )
            @test length(spl) == 14
            push!(psi_vals, spl[6])
            push!(conf_vals, spl[7])
         end
         println(stderr,"")


         out = IOBuffer(read=true,write=true)
         bs  = BufferedOutputStream(out)
         output_psi_header( bs )
         for g in 1:length(lib.graphs)
            name = lib.names[g]
            chr  = lib.info[g].name
            strand = lib.info[g].strand ? '+' : '-'
            _process_events( bs, lib.graphs[g], bam_quant.quant[g], (name,chr,strand), readlen=20 )
         end
         flush(bs)
         seek(out,0)
         i = 1
         for l in eachline(out)
            println(stderr,l)
            spl = split( l, '\t' )
            @test spl[6] == psi_vals[i]
            @test spl[7] == conf_vals[i]
            i += 1
         end

         # Verify TS tandem-start path-to-node ordering (node n must match path starting at n)
         ts_sg = SpliceGraph{2}(
            CoordInt[1, 101, 201],
            CoordInt[100, 200, 300],
            CoordInt[100, 100, 100],
            EdgeType[EDGETYPE_SL, EDGETYPE_SL, EDGETYPE_SL, EDGETYPE_RS],
            DNAMer{2}[DNAMer("AA") for _ in 1:4],
            DNAMer{2}[DNAMer("AA") for _ in 1:4],
            BitSet[BitSet([1, 2, 3]), BitSet([2, 3]), BitSet([3])],
            String["tx1", "tx2", "tx3"],
            ExpectedGC[zeros(Float16, 21) for _ in 1:3],
            LongSequence{DNAAlphabet{4}}("A"^300),
         )
         ts_quant = SpliceGraphQuant{SGAlignSingle, DefaultCounter}(ts_sg)
         ts_quant.leng .= 100.0
         # All reads only in downstream node 3 => downstream TSS (node 3) must get PSI = 1.0,
         # upstream TSS (node 1) must get PSI = 0.0, and each row's Inc_Paths must start at its node!
         ts_quant.node[3] = DefaultCounter(50.0)
         ts_psi, ts_utr, ts_ambig, ts_len = process_tandem_utr(ts_sg, ts_quant, NodeInt(1), TXST_MOTIF, 50)
         @test !isnull(ts_psi)
         ts_buf = IOBuffer(read=true, write=true)
         ts_bs = BufferedOutputStream(ts_buf)
         output_utr(ts_bs, map(x -> round(x, digits=4), get(ts_psi)), ts_utr, sum(ts_utr) + sum(ts_ambig), "TS", ts_sg, 1, ("TS_GENE", "chr1", '+'))
         flush(ts_bs)
         seek(ts_buf, 0)
         ts_lines = collect(eachline(ts_buf))
         @test length(ts_lines) == 3
         row1 = split(ts_lines[1], '\t')
         row2 = split(ts_lines[2], '\t')
         row3 = split(ts_lines[3], '\t')
         @test row1[2] == "1" && parse(Float64, row1[6]) == 0.0 && startswith(row1[12], "1-2-3:")
         @test row2[2] == "2" && parse(Float64, row2[6]) == 0.0 && startswith(row2[12], "2-3:")
         @test row3[2] == "3" && parse(Float64, row3[6]) == 1.0 && startswith(row3[12], "3:")
      end

      @testset "Diff and Delta" begin
         # 1. Beta-Binomial MLE & LRT on strong 2-group difference (including fractional counts)
         k = Float64[80.5, 76.2, 82.0, 20.1, 18.4, 21.5]
         n = Float64[100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
         r = n .- k
         X = [1.0 1.0; 1.0 1.0; 1.0 1.0; 1.0 0.0; 1.0 0.0; 1.0 0.0]
         beta, ll_full, H, conv = fit_bb_beta(k, r, n, X, 50.0)
         @test conv
         @test isapprox(_logistic(beta[1]), 0.20, atol=0.02)
         @test isapprox(_logistic(beta[1] + beta[2]), 0.80, atol=0.02)
         _, ll_null, _, _ = fit_bb_beta(k, r, n, X[:, 1:1], 50.0)
         @test 2.0 * (ll_full - ll_null) > 30.0

         # 2. Complete separation (1.0 vs 0.0 PSI) converges cleanly without NaN/Inf
         k_sep = Float64[50.0, 40.0, 60.0, 0.0, 0.0, 0.0]
         n_sep = Float64[50.0, 40.0, 60.0, 50.0, 40.0, 60.0]
         r_sep = n_sep .- k_sep
         b_sep, ll_sep, _, _ = fit_bb_beta(k_sep, r_sep, n_sep, X, 20.0)
         @test all(isfinite, b_sep) && isfinite(ll_sep)
         @test _logistic(b_sep[1] + b_sep[2]) > 0.98
         @test _logistic(b_sep[1]) < 0.02

         # 3. Benjamini-Hochberg FDR (including NaN preservation)
         fdr = bh_fdr([0.001, 0.04, NaN, 0.02, 0.5])
         @test isapprox(fdr[1], 0.004, atol=1e-12)
         @test isapprox(fdr[4], 0.04, atol=1e-12)
         @test isapprox(fdr[2], 0.053333333333, atol=1e-9)
         @test isnan(fdr[3])
         @test isapprox(fdr[5], 0.5, atol=1e-12)

         # 4. Formula parsing (additive, interaction, * expansion)
         @test parse_formula_terms("~ batch + donor + condition") == [["batch"], ["donor"], ["condition"]]
         @test parse_formula_terms("~ batch * condition") == [["batch"], ["condition"], ["batch", "condition"]]

         # 5. End-to-end covariate-adjusted Beta-Binomial GLM on synthetic batch-confounded .psi.gz files
         #    6 samples across 2 batches (B1, B2) and 2 conditions (ctrl, trt), plus continuous covariate `age`:
         #    - Event 1 (strong true condition effect + batch shift):
         #        B1 ctrl ≈ 0.20, B1 trt ≈ 0.70; B2 ctrl ≈ 0.40, B2 trt ≈ 0.85
         #    - Event 2 (pure batch effect, ZERO condition effect):
         #        B1 ctrl ≈ 0.20, B1 trt ≈ 0.20; B2 ctrl ≈ 0.80, B2 trt ≈ 0.80
         #    - Events 3..24 (null background for EB dispersion prior estimation)
         tmpdir = mktempdir()
         psi_files = String[]
         meta_rows = [
            ("s1", "s1.psi.gz", "ctrl", "B1", 30.0),
            ("s2", "s2.psi.gz", "ctrl", "B1", 35.0),
            ("s3", "s3.psi.gz", "ctrl", "B2", 50.0),
            ("s4", "s4.psi.gz", "trt",  "B1", 32.0),
            ("s5", "s5.psi.gz", "trt",  "B2", 52.0),
            ("s6", "s6.psi.gz", "trt",  "B2", 58.0),
         ]
         ev1_psi = [0.20, 0.22, 0.40, 0.72, 0.84, 0.86]
         ev2_psi = [0.20, 0.21, 0.80, 0.19, 0.79, 0.81]
         rng_glm = MersenneTwister(2026)
         bg_psi = [clamp.(0.5 .+ 0.05 .* randn(rng_glm, 6), 0.05, 0.95) for _ in 1:22]

         for (si, (sname, fname, cond, batch, age)) in enumerate(meta_rows)
            fpath = joinpath(tmpdir, fname)
            push!(psi_files, fpath)
            io_p = open(fpath, "w")
            st_p = ZlibDeflateOutputStream(io_p)
            write(st_p, "Gene\tNode\tCoord\tStrand\tType\tPsi\tCI_Width\tCI_Lo,Hi\tTotal_Reads\tComplexity\tEntropy\tInc_Paths\tExc_Paths\tEdges\n")
            write(st_p, "GENE1\t2\tchr1:100-200\t+\tCE\t$(ev1_psi[si])\t0.1\t0.1,0.2\t100.0\tK1\t0.8\t1\t1\t1-2:50\n")
            write(st_p, "GENE2\t2\tchr1:300-400\t+\tCE\t$(ev2_psi[si])\t0.1\t0.1,0.2\t100.0\tK1\t0.8\t1\t1\t1-2:50\n")
            for bi in 1:length(bg_psi)
               write(st_p, "BG$(bi)\t2\tchr1:$(500+bi*100)-$(550+bi*100)\t+\tCE\t$(round(bg_psi[bi][si], digits=4))\t0.1\t0.4,0.6\t80.0\tK1\t0.9\t1\t1\t1-2:40\n")
            end
            close(st_p); close(io_p)
         end

         design_tsv = joinpath(tmpdir, "design.tsv")
         open(design_tsv, "w") do io_d
            write(io_d, "sample\tfile\tcondition\tbatch\tage\n")
            for (sname, fname, cond, batch, age) in meta_rows
               write(io_d, "$sname\t$fname\t$cond\t$batch\t$age\n")
            end
         end

         design_cov = build_design_from_table(design_tsv, "~ batch + condition";
                                              test_spec="condition:trt_vs_ctrl", dir=tmpdir)
         @test design_cov.is_two_group
         @test size(design_cov.X) == (6, 3)

         out_cov = joinpath(tmpdir, "out_cov.diff.gz")
         evs_c, res_c, fdr_c = process_psi_files_glm(out_cov, design_cov; min_samp=2, min_reads=5)
         @test length(evs_c) == 24
         # Event 1: true condition effect detected with positive DeltaPsi (trt > ctrl) and tiny FDR
         @test res_c[1].tested
         @test res_c[1].delta_psi > 0.40
         @test fdr_c[1] < 1e-3
         # Event 2: pure batch effect completely controlled for by `~ batch + condition`!
         @test res_c[2].tested
         @test abs(res_c[2].delta_psi) < 0.05
         @test res_c[2].pvalue > 0.30

         # 6. Classic 2-group design (`build_two_group_design`) sign convention: DeltaPsi = Psi_A - Psi_B
         design_2g = build_two_group_design(psi_files[4:6], psi_files[1:3])
         out_2g = joinpath(tmpdir, "out_2g.diff.gz")
         _, res_2g, _ = process_psi_files_glm(out_2g, design_2g; min_samp=2, min_reads=5)
         @test res_2g[1].delta_psi > 0.50
         @test isapprox(res_2g[1].delta_psi, res_2g[1].psi_a - res_2g[1].psi_b, atol=1e-9)

         # 7. Continuous covariate test (`~ age`)
         design_age = build_design_from_table(design_tsv, "~ age"; test_spec="age", dir=tmpdir)
         @test design_age.is_continuous_test
         out_age = joinpath(tmpdir, "out_age.diff.gz")
         _, res_age, _ = process_psi_files_glm(out_age, design_age; min_samp=2, min_reads=5)
         @test res_age[2].tested && res_age[2].pvalue < 0.01 # Event 2 correlates strongly with B2 (older age)

         # 8. CLI end-to-end check on both new GLM (--design and -a/-b) and --legacy modes
         run(`$(Base.julia_cmd()) --project=.. ../bin/whippet-delta.jl --design $design_tsv --formula "~ batch + condition" --test condition:trt_vs_ctrl -d $tmpdir -o $(joinpath(tmpdir, "cli_cov"))`)
         @test isfile(joinpath(tmpdir, "cli_cov.diff.gz"))
         run(`$(Base.julia_cmd()) --project=.. ../bin/whippet-delta.jl -a $(join(psi_files[4:6], ',')) -b $(join(psi_files[1:3], ',')) --legacy -o $(joinpath(tmpdir, "cli_leg"))`)
         @test isfile(joinpath(tmpdir, "cli_leg.diff.gz"))
         rm(tmpdir, recursive=true, force=true)
      end
   end

   @testset "Executable testing" begin
      check_samtools() || error("Samtools is not properly installed, but is a dependency of BAM input!")
      run(`pwd`)
      run(`$(Base.julia_cmd()) ../bin/whippet-index.jl -h`)
      run(`$(Base.julia_cmd()) ../bin/whippet-quant.jl -h`)
      run(`$(Base.julia_cmd()) ../bin/whippet-pivot.jl -h`)
      run(`$(Base.julia_cmd()) ../bin/whippet-motif.jl -h`)
      run(`$(Base.julia_cmd()) ../bin/whippet-delta.jl -h`)
      run(`samtools view -b test_out.sam -o test_out.bam`)
      run(`samtools sort test_out.bam -o test_out.sort.bam`)
      run(`samtools index test_out.sort.bam`)
      run(`$(Base.julia_cmd()) ../bin/whippet-quant.jl test_out.sort.bam -x test_index.jls`)
      run(`$(Base.julia_cmd()) ../bin/whippet-quant.jl aber_inp.bam -x test_index.jls`)
      run(`$(Base.julia_cmd()) ../bin/whippet-quant.jl test.fq.gz -x test_index.jls`)
      run(`$(Base.julia_cmd()) ../bin/whippet-quant.jl test_real.fq.gz -x test_index.jls`)
       ## CHECK test.sort.bam, make 
   end
end
