using System;
using System.IO;
using System.Collections.Generic;
using NUnit;
using NUnit.Framework;
using Bio.BWA.MEM;
using Bio;
using Bio.Extensions;
using Bio.BWA;
using System.Linq;

namespace TestBWA
{
    [TestFixture]
    public static class TestBWA
    {
        const string fasta_name = "../../../../TestData/MT.fasta";

        [Test]
        [Category("BWA")]
        public static void TestFastaCreation ()
        {
            var exists = File.Exists (fasta_name);
            Assert.IsTrue (exists);
            using (BWA bwa = new BWA (fasta_name)) {
                Assert.Pass ();
            }                
        }

        [Test]
        [Category("BWA")]
        public static void TestForwardAlignment() {
            var seq = new Sequence (DnaAlphabet.Instance, "TCTACAAACCACAAAGACATTGGAACACTATACCTATTA");
            using (BWA bwa = new BWA (fasta_name)) {
                var aln = bwa.AlignSequence (seq);
                Assert.AreEqual (seq.Count.ToString () + "M", aln.CIGAR);
                Assert.AreEqual (5927, aln.Pos);
                Assert.AreEqual ("MT", aln.RName);
            }
        }

        [Test]
        [Category("BWA")]
        public static void TestReverseAlignment() {
            var seqstr = "TAATAGGTATAGTGTTCCAATGTCTTTGTGGTTTGTAGA";
            var seq = new Sequence (DnaAlphabet.Instance, seqstr);
            using (BWA bwa = new BWA (fasta_name)) {
                var aln = bwa.AlignSequence (seq);
                Assert.AreEqual (seq.Count.ToString () + "M", aln.CIGAR);
                Assert.AreEqual (5927, aln.Pos);
                Assert.AreEqual (seqstr.Length, aln.QuerySequence.Count);
                Assert.AreEqual ("TCTACAAACCACAAAGACATTGGAACACTATACCTATTA", aln.QuerySequence.ConvertToString ());
                Assert.AreEqual (16, (int)aln.Flag);
                Assert.AreEqual (39, aln.Metadata ["AS"]);
                Assert.AreEqual (0, aln.Metadata ["XS"]);
                Assert.AreEqual (0, aln.Metadata ["NM"]);
                Assert.AreEqual ("MT", aln.RName);
                Assert.IsTrue (aln.QuerySequence.IsMarkedAsReverseComplement ());
            }
        }
        [Test]
        [Category("BWA")]
        public static void TestSNPAlignment() {
            var seq = new Sequence (DnaAlphabet.Instance, "TCTACAAACCACAAAGACATTGGAAACTATACCTATTA");
            using (BWA bwa = new BWA (fasta_name)) {
                var aln = bwa.AlignSequence (seq);
                Assert.AreEqual ("25M1D13M", aln.CIGAR);
                Assert.AreEqual (5927, aln.Pos);
                Assert.AreEqual ("MT", aln.RName);
                Assert.AreEqual (1, aln.Metadata ["NM"]);
            }
        }

        [Test]
        [Category("BWA")]
        public static void TestClippedReverseAlignment() {
            var seq = new Sequence (DnaAlphabet.Instance, "GGGGGGTAATAGGTATAGTGTTCCAATGTCTTTGTGGTTTGTAGAGGGGGGGG");
            using (BWA bwa = new BWA (fasta_name)) {
                var aln = bwa.AlignSequence (seq);
                Assert.AreEqual ("7S40M6S", aln.CIGAR);
                Assert.AreEqual ("CCCCCCCCTCTACAAACCACAAAGACATTGGAACACTATACCTATTACCCCCC", aln.QuerySequence.ConvertToString ());
                Assert.AreEqual (5926, aln.Pos);
                Assert.AreEqual ("MT", aln.RName);
                Assert.AreEqual (Bio.IO.SAM.SAMFlags.QueryOnReverseStrand, aln.Flag);
                Assert.IsTrue (aln.QuerySequence.IsMarkedAsReverseComplement ());
            }
        }

        [Test]
        [Category("BWA")]
        public static void TestClippedPairwiseAlignment() {
            var seq = new Sequence (DnaAlphabet.Instance, "GGGGGGTAATAGGTATAGTGTTCCAATGTCTTTGTGGTTTGTAGAGGGGGGGG");

            using (var bwa = new BWAPairwiseAligner (fasta_name, false)) {
                var aln = bwa.AlignRead (seq);
                Assert.AreEqual (aln.PairwiseAlignedSequences.Count, 1);
                var refseq = aln.PairwiseAlignedSequences[0].FirstSequence;
                var queryseq = aln.PairwiseAlignedSequences [0].SecondSequence;
                Assert.AreEqual ("CTCTACAAACCACAAAGACATTGGAACACTATACCTATTA", queryseq.ConvertToString ());
                Assert.IsTrue (queryseq.IsMarkedAsReverseComplement ());
                Assert.AreEqual ("CTCTACAAACCACAAAGACATTGGAACACTATACCTATTA", refseq.ConvertToString ());
                var baln = aln as BWAPairwiseAlignment;
                Assert.AreEqual (5926, baln.AlignedSAMSequence.Pos);
                Assert.AreEqual(queryseq, baln.AlignedQuerySeq);
                Assert.AreEqual(refseq, baln.AlignedRefSeq);
            }
        }

        [Test]
        [Category("BWA")]
        public static void TestDeletionAlignment() {
            var seq = new Sequence (DnaAlphabet.Instance, "TCTACAAACCACAAAGACATTGGAAACTATACCTATTA");
            using (var bwa = new BWAPairwiseAligner (fasta_name, false)) {
                var aln = bwa.AlignRead (seq);
                var baln = aln as BWAPairwiseAlignment;
                Assert.AreEqual ("25M1D13M", baln.AlignedSAMSequence.CIGAR);
                Assert.AreEqual (aln.PairwiseAlignedSequences.Count, 1);
                var refseq = aln.PairwiseAlignedSequences[0].FirstSequence;
                var queryseq = aln.PairwiseAlignedSequences [0].SecondSequence;
                Assert.AreEqual ('-', queryseq [25]);
                Assert.AreEqual ("TCTACAAACCACAAAGACATTGGAA-ACTATACCTATTA", queryseq.ConvertToString ());
                Assert.AreEqual ("TCTACAAACCACAAAGACATTGGAACACTATACCTATTA", refseq.ConvertToString ());
                Assert.AreEqual (5927, baln.AlignedSAMSequence.Pos);
                Assert.AreEqual(queryseq, baln.AlignedQuerySeq);
                Assert.AreEqual(refseq, baln.AlignedRefSeq);
            }            
        }

        // Can revitalize this later, is for a secondary alignment, so doesn't apply when we only return the best
        #if FALSE
        [Test]
        [Category("BWA")]
        public static void TestCCSReads() {
            var seq =  "GAAGCTACTAGTCCTCAGCAAGCTTGTGCGTCGCTCAAAAAGCTGCGCTCGAAAAAAAAAAGTCGTCTGTCTAGATGATGTGCCCCCCCCCCGTATATGTATCCCCCAGTGTATGAGCATTCTAGAGGATCCCCGGGTCTCTCTCAAGAATGCTCATACACTGGGGGATACATATACGGGGGGGGGGCACATCATCTAGACAGACGACTTTTTTTTTCGAGCGCAGCTTTTTGAGCGACGCACAAGCTTGCTGAGGACTAGTAGCTTC";
            var qual = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~3~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~E~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~J~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
            var seqq = new QualitativeSequence (DnaAlphabet.Instance, FastQFormatType.Sanger, 
                           seq.Select (x => (byte)x).ToArray (),
                           qual.Select (p => (int)p).ToArray (),
                false);
            var fname = "../../../../TestData/References.fna";
            using (var bwa = new BWAPairwiseAligner (fname, false)) {
                var aln = bwa.AlignRead (seqq);
                var baln = aln as BWAPairwiseAlignment;
                Assert.AreEqual (14, baln.AlignedSAMSequence.Pos);
                Assert.AreEqual ("146S62M1D60M", baln.AlignedSAMSequence.CIGAR);
                Assert.AreEqual (aln.PairwiseAlignedSequences.Count, 1);
                var refseq = aln.PairwiseAlignedSequences[0].FirstSequence;
                var queryseq = aln.PairwiseAlignedSequences [0].SecondSequence;
                Assert.AreEqual ('-', queryseq [62]);
                Assert.AreEqual ("AGAATGCTCATACACTGGGGGATACATATACGGGGGGGGGGCACATCATCTAGACAGACGAC-TTTTTTTTTCGAGCGCAGCTTTTTGAGCGACGCACAAGCTTGCTGAGGACTAGTAGCTTC", 
                                    queryseq.ConvertToString ());
                Assert.AreEqual ("AGAATGCTCATACACTGGGGGATACATATACGGGGGGGGGGCACATCATCTAGACAGACGACTTTTTTTTTTCGAGCGCAGCTTTTTGAGCGACGCACAAGCTTGCTGAGGACTAGTAGCTTC", 
                    refseq.ConvertToString ());
                Assert.AreEqual(queryseq, baln.AlignedQuerySeq);
                Assert.AreEqual(refseq, baln.AlignedRefSeq);
            }   
        }
        #endif

        [Test]
        [Category("BWA")]
        public static void TestCCSReads() {
            var seq =  "GAAGCTACTAGTCCTCAGCAAGCTTGTGCGTCGCTCAAAAAGCTGCGCTCGAAAAAAAAAAGTCGTCTGTCTAGATGATGTGCCCCCCCCCCGTATATGTATCCCCCAGTGTATGAGCATTCTAGAGGATCCCCGGGTCTCTCTCAAGAATGCTCATACACTGGGGGATACATATACGGGGGGGGGGCACATCATCTAGACAGACGACTTTTTTTTTCGAGCGCAGCTTTTTGAGCGACGCACAAGCTTGCTGAGGACTAGTAGCTTC";
            var qual = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~3~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~E~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~J~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
            var seqq = new QualitativeSequence (DnaAlphabet.Instance, FastQFormatType.Sanger, 
                seq.Select (x => (byte)x).ToArray (),
                qual.Select (p => (int)p).ToArray (),
                false);
            var fname = "../../../../TestData/References.fna";
            using (var bwa = new BWAPairwiseAligner (fname, false)) {
                var aln = bwa.AlignRead (seqq);
                var baln = aln as BWAPairwiseAlignment;
                Assert.AreEqual (1, baln.AlignedSAMSequence.Pos);
                Assert.AreEqual ("132S136M", baln.AlignedSAMSequence.CIGAR);
                Assert.AreEqual (aln.PairwiseAlignedSequences.Count, 1);
                var refseq = aln.PairwiseAlignedSequences[0].FirstSequence;
                var queryseq = aln.PairwiseAlignedSequences [0].SecondSequence;
                var exp_seq = "CCGGGGATCCTCTAGAATGCTCATACACTGGGGGATACATATACGGGGGGGGGGCACATCATCTAGACAGACGACTTTTTTTTTTCGAGCGCAGCTTTTTGAGCGACGCACAAGCTTGCTGAGGACTAGTAGCTTC";
                Assert.AreEqual (exp_seq, 
                    queryseq.ConvertToString ());
                Assert.AreEqual (exp_seq, 
                    refseq.ConvertToString ());
                Assert.AreEqual(queryseq, baln.AlignedQuerySeq);
                Assert.AreEqual(refseq, baln.AlignedRefSeq);
            }   
        }

        [Test]
        [Category("BWA")]
        public static void Test6S1I105M1I32MCigar() {
            /* This read comes back with a cigar having both an insertion and a match, this combines them here so 6S1I -> 7S */
            var seq =  "GAAGCTACTAGTCCTCAGCAAGCTTGTGCGTCCGCTCAAAAAGCTGCGCTCGAAAAAAAAAAGTCGTCTGTCTAGATGATGTGCCCCCCCCCCGTATATGTATCCCCCAGTGTATGAGCATTCTAGAGGATCCCCGGGATGCTCT";
            var qual = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~3~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
            var seqq = new QualitativeSequence (DnaAlphabet.Instance, FastQFormatType.Sanger, 
                seq.Select (x => (byte)x).ToArray (),
                qual.Select (p => (int)p).ToArray (),
                false);
            var fname = "../../../../TestData/References.fna";
            using (var bwa = new BWAPairwiseAligner (fname, false)) {
                var aln = bwa.AlignRead (seqq);
                var baln = aln as BWAPairwiseAlignment;
                Assert.AreEqual (0, baln.AlignedSAMSequence.Pos);
                // Currently fixed when generating the alignment, but not the read.
                Assert.AreEqual ("6S1I105M1I32M", baln.AlignedSAMSequence.CIGAR);
                Assert.AreEqual (aln.PairwiseAlignedSequences.Count, 1);
                Assert.AreEqual (Bio.IO.SAM.SAMFlags.QueryOnReverseStrand, baln.AlignedSAMSequence.Flag);
                var refseq = aln.PairwiseAlignedSequences[0].FirstSequence;
                var queryseq = aln.PairwiseAlignedSequences [0].SecondSequence;
                var exp_seq = "CCCGGGGATCCTCTAGAATGCTCATACACTGGGGGATACATATACGGGGGGGGGGCACATCATCTAGACAGACGACTTTTTTTTTTCGAGCGCAGCTTTTTGAGCGGACGCACAAGCTTGCTGAGGACTAGTAGCTTC";
                Assert.AreEqual (exp_seq, 
                    queryseq.ConvertToString ());
                exp_seq = "CCCGGGGATCCTCTAGAATGCTCATACACTGGGGGATACATATACGGGGGGGGGGCACATCATCTAGACAGACGACTTTTTTTTTTCGAGCGCAGCTTTTTGAGC-GACGCACAAGCTTGCTGAGGACTAGTAGCTTC";

                Assert.AreEqual (exp_seq, 
                    refseq.ConvertToString ());
                Assert.AreEqual(queryseq, baln.AlignedQuerySeq);
                Assert.AreEqual(refseq, baln.AlignedRefSeq);
            }   
        }
    }
}

