using System;
using System.IO;
using System.Collections.Generic;
using NUnit;
using NUnit.Framework;
using Bio.BWA.MEM;
using Bio;

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
            var seq = new Sequence (DnaAlphabet.Instance, "TAATAGGTATAGTGTTCCAATGTCTTTGTGGTTTGTAGA");
            using (BWA bwa = new BWA (fasta_name)) {
                var aln = bwa.AlignSequence (seq);
                Assert.AreEqual (seq.Count.ToString () + "M", aln.CIGAR);
                Assert.AreEqual (5927, aln.Pos);
                Assert.AreEqual ((int)aln.Flag, 250);
                //Assert.AreEqual(aln.
                Assert.AreEqual ("MT", aln.RName);
            }

        }
    }
}

