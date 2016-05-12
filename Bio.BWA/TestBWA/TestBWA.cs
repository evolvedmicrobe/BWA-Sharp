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
        static string base_dir = "../../../..";
        //static string base_dir = "/Users/nigel/git/BWA-Sharp";
        static string fasta_name = Path.Combine(base_dir, "TestData/MT.fasta");
        static string screen_refs = Path.Combine(base_dir, "TestData/References.fna");

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
           
            using (var bwa = new BWAPairwiseAligner (screen_refs, false)) {
                var aln = bwa.AlignRead (seqq);
                var baln = aln as BWAPairwiseAlignment;
                Assert.AreEqual (0, baln.AlignedSAMSequence.Pos);
                Assert.AreEqual ("131S137M", baln.AlignedSAMSequence.CIGAR);
                Assert.AreEqual (aln.PairwiseAlignedSequences.Count, 1);
                var refseq = aln.PairwiseAlignedSequences[0].FirstSequence;
                var queryseq = aln.PairwiseAlignedSequences [0].SecondSequence;
                var exp_seq = "CCCGGGGATCCTCTAGAATGCTCATACACTGGGGGATACATATACGGGGGGGGGGCACATCATCTAGACAGACGACTTTTTTTTTTCGAGCGCAGCTTTTTGAGCGACGCACAAGCTTGCTGAGGACTAGTAGCTTC";
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
        public static void TestQVsAreFaithful() {
            var seq =  "GAAGCTACTAGTCCTCAGCAAGCTTGTGCGTCGCTCAAAAAGCTGCGCTCGAAAAAAAAAAGTCGTCTGTCTAGATGATGTGCCCCCCCCCCGTATATGTATCCCCCAGTGTATGAGCATTCTAGAGGATCCCCGGGTCTCTCTCAAGAATGCTCATACACTGGGGGATACATATACGGGGGGGGGGCACATCATCTAGACAGACGACTTTTTTTTTCGAGCGCAGCTTTTTGAGCGACGCACAAGCTTGCTGAGGACTAGTAGCTTC";
            var qual = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~3~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~E~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~J~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
            var seqq = new QualitativeSequence (DnaAlphabet.Instance, FastQFormatType.Sanger, 
                seq.Select (x => (byte)x).ToArray (),
                qual.Select (p => (int)p).ToArray (),
                false);
            using (var bwa = new BWAPairwiseAligner (screen_refs, false)) {
                var aln = bwa.AlignRead (seqq);
                var baln = aln as BWAPairwiseAlignment;
                Assert.AreEqual (0, baln.AlignedSAMSequence.Pos);
                Assert.AreEqual ("131S137M", baln.AlignedSAMSequence.CIGAR);
                Assert.AreEqual (aln.PairwiseAlignedSequences.Count, 1);
                var refseq = aln.PairwiseAlignedSequences[0].FirstSequence;
                var queryseq = aln.PairwiseAlignedSequences [0].SecondSequence;
                var exp_seq = "CCCGGGGATCCTCTAGAATGCTCATACACTGGGGGATACATATACGGGGGGGGGGCACATCATCTAGACAGACGACTTTTTTTTTTCGAGCGCAGCTTTTTGAGCGACGCACAAGCTTGCTGAGGACTAGTAGCTTC";
                Assert.AreEqual (exp_seq, 
                    queryseq.ConvertToString ());
                Assert.AreEqual (exp_seq, 
                    refseq.ConvertToString ());
                Assert.AreEqual(queryseq, baln.AlignedQuerySeq);
                Assert.AreEqual(refseq, baln.AlignedRefSeq);
            }   
        }

        //This was broken in the old BWA MEM
        [Test]
        [Category("BWA")]
        public static void TestCCSReadWithMismatchAtEnd() {
            var seq =  "CTGGCAACTAATTCAGTCCAGTAAATATCCTCAATAGGGAATAATATATGCTTTCCATTCCATCGGGAAAAAGTTTTGTTCAACACACCAAGCTCAATCAACTCACTAATGTATGGGAATTGTTTTGATGTAACCACATACTTCCTGCCTTCATTAAGGGCTGCGCACAAAACCATAGATTGCTCTTCTGTAAGGTTTTGAATTACTGATCGCACTTTATCGTTTTGCATCTTAATGCGTTTTCTTAGCTTAAATCGCTTATATCTGGCGCTGGCAATAGCTGATAATCGATGCACATTAATTGCTAGCGAAAATGCAAGAGCAAAGACGAAAACATGCCACACATGAGGAATACCGATTCTCTCATTAACATATTCAGGCCAGTTATCTGGGCTTAAAAGCAGAAGTCCAACCCAGATAACGATCATATACATGGTTCTCTCCAGAGGTTCATTACTGAACACTCGTCCGAGAATAACGAGTGGATCCATTTCTATACTCATCAAACTGTAGGGGTTGTAATAGTTTATCCGATTTCTCGCTGTAGGGGTACACGAGAACCACCGAGCCTGATGTGGTTAAAAGACAGGCACAATCTTTACTACCGCAATCCACTATTTAAGGTGATATATGGAAGAAG";
            var qual = "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK@KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK";
            var seqq = new QualitativeSequence (DnaAlphabet.Instance, FastQFormatType.Sanger, 
                seq.Select (x => (byte)x).ToArray (),
                qual.Select (p => (int)p).ToArray (),
                false);

            using (var bwa = new BWA (screen_refs)) {
                var aln = bwa.AlignSequence (seqq);
                Assert.AreEqual (34347, aln.Pos);
                Assert.AreEqual ("640M", aln.CIGAR);
                Assert.AreEqual(0, aln.Metadata["NM"]);
            }

            using (var bwa = new BWAPairwiseAligner (screen_refs)) {
                var aln = bwa.AlignRead (seqq);
                var baln = aln as BWAPairwiseAlignment;
                // Assert.AreEqual (34348, baln.AlignedSAMSequence.Pos);
                Assert.AreEqual ("640M", baln.AlignedSAMSequence.CIGAR);
                Assert.AreEqual (aln.PairwiseAlignedSequences.Count, 1);
                var refseq = aln.PairwiseAlignedSequences[0].FirstSequence;
                var queryseq = aln.PairwiseAlignedSequences [0].SecondSequence;
                var exp_seq = "CTTCTTCCATATATCACCTTAAATAGTGGATTGCGGTAGTAAAGATTGTGCCTGTCTTTTAACCACATCAGGCTCGGTGGTTCTCGTGTACCCCTACAGCGAGAAATCGGATAAACTATTACAACCCCTACAGTTTGATGAGTATAGAAATGGATCCACTCGTTATTCTCGGACGAGTGTTCAGTAATGAACCTCTGGAGAGAACCATGTATATGATCGTTATCTGGGTTGGACTTCTGCTTTTAAGCCCAGATAACTGGCCTGAATATGTTAATGAGAGAATCGGTATTCCTCATGTGTGGCATGTTTTCGTCTTTGCTCTTGCATTTTCGCTAGCAATTAATGTGCATCGATTATCAGCTATTGCCAGCGCCAGATATAAGCGATTTAAGCTAAGAAAACGCATTAAGATGCAAAACGATAAAGTGCGATCAGTAATTCAAAACCTTACAGAAGAGCAATCTATGGTTTTGTGCGCAGCCCTTAATGAAGGCAGGAAGTATGTGGTTACATCAAAACAATTCCCATACATTAGTGAGTTGATTGAGCTTGGTGTGTTGAACAAAACTTTTTCCCGATGGAATGGAAAGCATATATTATTCCCTATTGAGGATATTTACTGGACTGAATTAGTTGCCAG";
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
            using (var bwa = new BWAPairwiseAligner (screen_refs, false)) {
                var aln = bwa.AlignRead (seqq);
                var baln = aln as BWAPairwiseAlignment;
                Assert.AreEqual (0, baln.AlignedSAMSequence.Pos);
                // Currently fixed when generating the alignment, but not the read.
                Assert.AreEqual ("7S105M1I32M", baln.AlignedSAMSequence.CIGAR);
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

        [Test]
        [Category("BWA")]
        public static void TestAddingNodesThatOverlapsWithTwoExistingNodes() {
            /* Here we test that a node that ovelaps with a region when added, does not generate 
             * three nodes, but only two */
            RegionTree rt = new RegionTree ();
            string chrm = "1";
            var left = new Region (chrm, 1, 5);
            var right = new Region (chrm, 10, 20);
            var middle = new Region (chrm, 3, 15);
            rt.Add (left);
            rt.Add (right);
            rt.Add (middle);
            Assert.AreEqual (1, rt.Count);

            var node = rt.InOrderTraversal ().ToList ().First();
            Assert.AreEqual (1, node.Start);
            Assert.AreEqual (20, node.End);
            Assert.AreEqual (chrm, node.Reference);


            // Now for a more complicated example
            var excludedRegion = new Region (chrm, 100, 200);
            rt = new RegionTree ();
            Random rv = new Random (56448638);
            for (int i = 0; i < 500; i++) {
                var s = rv.Next (1, 1000);
                var e = s + rv.Next (1, 5);
                var nr = new Region (chrm, s, e);
                if (nr.CompareTo (excludedRegion) != 0) {
                    rt.Add (nr);
                }
            }
            var allnodes = rt.InOrderTraversal ().ToList ();
            var pairs = Enumerable.Zip (allnodes.Take (allnodes.Count - 1), allnodes.Skip (1), (arg1, arg2) => new Tuple<Region, Region>(arg1, arg2)).ToList();
            // Verify they return in sorted order
            var sortedAndNoOverlap = !pairs.Any (x => x.Item1.CompareTo (x.Item2) >=0);
            Assert.True (sortedAndNoOverlap);
            // Does indexing count match actual count?
            Assert.AreEqual (allnodes.Count, rt.Count);

            // Now add one
            rt.Add (new Region(chrm, excludedRegion.Start - 20, excludedRegion.End + 20));
            Assert.AreEqual (156, rt.Count);

        }
    }
}

