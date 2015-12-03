using System;
using System.Linq;
using System.IO;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;

using Bio;
using Bio.Algorithms.Alignment;
using Bio.IO.FastA;
using Bio.IO.SAM;
using Bio.Variant;
using Bio.BWA.MEM;

namespace Bio.BWA
{
    /// <summary>
    /// A class to create pairwise alignments using BWA to do the alignment.  The BWA class returns SAM records from an alignment, while the 
    /// BWAPairwiseAligner goes through the additional step of adding in the reference genome and creating an IPairwiseAlignment object.
    /// 
    /// Useful for cases where we want more information than just what the SAM record provides.
    /// </summary>
    public class BWAPairwiseAligner : IDisposable
    {
        /// <summary>
        /// Determines the breakpoint between using a memory mapped file and holding stuff in RAM
        /// </summary>
        const int MAX_CACHED_FASTA_SIZE = 10000000;
        Bio.BWA.MEM.BWA bwa; 
        FastaSubSequenceProvider fi;
        RegionTree tree;
        int UnAlignedReads =0;
        int AlignedReads = 0;
        string fasta_name;
        bool recordCoverageIntervals;

        public BWAPairwiseAligner (string fasta_name, bool recordCoverageIntervals=false, bool usePacBioOptions=false)
        {
            this.fasta_name = fasta_name;
            this.recordCoverageIntervals = recordCoverageIntervals;
            if (recordCoverageIntervals) {
                tree = new RegionTree ();
            }
            bwa = new Bio.BWA.MEM.BWA (fasta_name, usePacBioOptions);

            var finfo = new FileInfo (fasta_name);
            if (finfo.Length > MAX_CACHED_FASTA_SIZE) {
                fi = new LargeFastaSubSequenceProvider (fasta_name);
            } else {
                fi = new SmallFastaSubSequenceProvider (fasta_name);
            }
        }

        public IPairwiseSequenceAlignment AlignRead(ISequence read)
        {
            var aln = bwa.AlignSequence (read);
            if (aln == null) {
                Interlocked.Increment (ref UnAlignedReads);
                return null;
            }
            else {
                Interlocked.Increment (ref AlignedReads);
                var refSeq = fi.GetReferenceSection (aln.RName, aln.Pos, aln.RefEndPos + 1);
                var expanded = expandAlignment (refSeq, aln);
#if OLD_BWA
                if (clipEnds) {
                    CorrectEndMismatches (ref expanded, ref aln);
                }
#endif
                var qseq = new QualitativeSequence(DnaAlphabet.Instance, FastQFormatType.Sanger, 
                    expanded.Item2.Select(x => x.BP).ToArray(),
                    expanded.Item2.Select(z=> (int) z.QV).ToArray(), false);
                var rseq = new Sequence (DnaAlphabet.Instance, expanded.Item1.ToArray (), false);

                var full_alignment = new BWAPairwiseAlignment (aln, rseq, qseq);
                                  
                if (recordCoverageIntervals) {
                    lock (tree) {
                        var reg = new Region (aln.RName, aln.Pos, aln.RefEndPos);
                        tree.Add (reg);
                    }
                }
                return full_alignment;
            }
        }

        Tuple<List<byte>, List<BPandQV>> expandAlignment(byte[] refSeq, SAMAlignedSequence aln)
        {
            if (string.IsNullOrWhiteSpace (aln.CIGAR) || aln.CIGAR.Equals ("*")) {
                throw new Exception ();
            }

            var elements = CigarUtils.GetCigarElements (aln.CIGAR);
            // This used to be necessary but was fixed in later BWA versions it seems.
            //elements = CigarUtils.FilterLeadingAndTrailingInsertions (elements);

            int length = CigarUtils.GetAlignmentLengthAfterClipping (elements);
            var q = new List<BPandQV>(length);
            var r = new List<byte>(length);

            int curRef = 0;
            int curQuery = 0; //location on query
            var qseq = aln.QuerySequence as QualitativeSequence;
            bool isQual = true;
            if (qseq == null) {
                isQual = false;
            }
            var qstr = aln.QuerySequence;

            // Account for soft clipping in offset
            foreach (var e in elements) {
                char ch = e.Operation;
                // Check below is currently also implemented by GetAlignmentLengthAfterClipping
                // can remove if profiling shows expensive
                if (ch == 'P' || ch == 'N' || ch=='H') {
                    throw new Exception ("Not built to handle clipping yet");
                }
                int len = e.Length;
                bool addRef = false;
                bool addQuery = false;
                if (ch == 'M' || ch == '=' || ch == 'X') {
                    addRef = true;
                    addQuery = true;
                } else if (ch == 'D') {
                    addRef = true;
                } else if (ch == 'I') {
                    addQuery = true;
                } else if (ch == 'S') {
                    curQuery += len;
                    continue;
                }
                for (int j = 0; j < len; j++) {
                    if (addRef) {
                        r.Add (refSeq [curRef]);
                        curRef++;
                    } else {
                        r.Add ((byte)'-');
                    }
                    if (addQuery) {
                        var bp = (byte)qstr [curQuery];
                        var qv = isQual ? (byte) qseq.GetQualityScore(curQuery) : (byte)0;
                        q.Add (new BPandQV (bp, qv));
                        curQuery++;
                    } else {
                        q.Add (new BPandQV ((byte)'-', 0));
                    }
                }

            }
            return new Tuple<List<byte>, List<BPandQV>> (r, q);
        }

        #if OLD_BWA
        /// <summary>
        /// BWA MEM will report bases with a mismatch at the end rather than soft or hard clipping them.
        /// Super annoying, but to maintain compatability with the default I fix this here.
        /// TODO: FIX THIS!
        /// </summary>
        /// <param name="expanded">Expanded.</param>
        /// <param name="aln">Aln.</param>
        void CorrectEndMismatches( ref Tuple<List<byte>, List<BPandQV>> expanded, ref SAMAlignedSequence aln) {
            var refseq = expanded.Item1;
            var qseq = expanded.Item2;
            // Fix mismatch in start
            if (refseq [0] != qseq [0].BP) {

                // Make sure we only have one mismatch
                if (refseq.Count > 1 && qseq.Count > 1 &&
                    refseq [1] != qseq [1].BP) {
                    throw new BWAException ("BWA returned two mismatches at the start of the alignment.  This is a bug that should be reported.");
                }

                // Fix the cigar
                var elements = CigarUtils.GetCigarElements (aln.CIGAR);
                if (elements.First ().Operation == CigarOperations.SOFT_CLIP &&
                    elements.Count > 2 &&
                    elements [1].Operation == CigarOperations.ALN_MATCH) {
                    var op = elements [0];
                    op.Length += 1;
                    elements [0] = op;
                    var op2 = elements [1];
                    op2.Length -= 1;
                    elements [1] = op2;
                } else if (elements.First ().Operation == CigarOperations.ALN_MATCH) {
                    var op = new CigarElement (CigarOperations.SOFT_CLIP, 1);
                    elements.Insert (0, op);
                    var op2 = elements [1];
                    op2.Length -= 1;
                    elements [1] = op2;
                } else {
                    throw new BWAException ("Found mismatch at start of alignment, but could not correct because cigar did not start with match or soft clip character");
                }

                SAMAlignedSequenceHeader header = new SAMAlignedSequenceHeader ();
                // Set invalid bin since it might no longer be accurate
                header.Pos = aln.Pos + 1;
                Console.WriteLine (header.Pos);
                header.CIGAR = CigarUtils.CreateCigarString (elements);
                header.Flag = aln.Flag;
                header.ISize = aln.ISize;
                header.MapQ = aln.MapQ;
                header.MPos = aln.MPos;
                header.MRNM = aln.MRNM;
                header.QName = aln.QName;
                header.RName = aln.RName;
                SAMAlignedSequence n_aln = new SAMAlignedSequence (header);
                n_aln.QuerySequence = aln.QuerySequence;
                aln = n_aln;
                var rnew = refseq.Skip (1).ToList ();
                var qnew = qseq.Skip (1).ToList ();
                expanded = new Tuple<List<byte>, List<BPandQV>> (rnew, qnew);
                refseq = expanded.Item1;
                qseq = expanded.Item2;

            } 
            // Fix end
            if (refseq.Last () != qseq.Last ().BP) {
                // Make sure we only have one mismatch
                if (refseq.Count > 1 && qseq.Count > 1 &&
                    refseq [refseq.Count - 2 ] != qseq [qseq.Count - 2].BP) {
                    throw new BWAException ("BWA returned two mismatches at the end of the alignment.  This is a bug that should be reported.");
                }

                // Fix the cigar
                var elements = CigarUtils.GetCigarElements (aln.CIGAR);
                if (elements.Last ().Operation == CigarOperations.SOFT_CLIP &&
                    elements.Count > 2 &&
                    elements [elements.Count - 2].Operation == CigarOperations.ALN_MATCH) {
                    var op = elements [elements.Count - 1];
                    op.Length += 1;
                    elements [elements.Count - 1] = op;
                    var op2 = elements [elements.Count - 2];
                    op2.Length -= 1;
                    elements [elements.Count - 2] = op2;
                } else if (elements.Last ().Operation == CigarOperations.ALN_MATCH) {
                    var op = new CigarElement (CigarOperations.SOFT_CLIP, 1);
                    elements.Add(op);
                    var op2 = elements [elements.Count - 2];
                    op2.Length -= 1;
                    elements [elements.Count - 2] = op2;
                } else {
                    throw new BWAException ("Found mismatch at start of alignment, but could not correct because cigar did not start with match or soft clip character");
                }

                SAMAlignedSequenceHeader header = new SAMAlignedSequenceHeader ();
                // Set invalid bin since it might no longer be accurate
                header.CIGAR = CigarUtils.CreateCigarString (elements);
                header.Flag = aln.Flag;
                header.ISize = aln.ISize;
                header.MapQ = aln.MapQ;
                header.MPos = aln.MPos;
                header.MRNM = aln.MRNM;
                header.QName = aln.QName;
                header.RName = aln.RName;
                SAMAlignedSequence n_aln = new SAMAlignedSequence (header);
                n_aln.QuerySequence = aln.QuerySequence;
                aln = n_aln;
                var rnew = refseq.Skip (1).ToList ();
                var qnew = qseq.Skip (1).ToList ();
                expanded = new Tuple<List<byte>, List<BPandQV>> (rnew, qnew);
            }
        }
        #endif


        public void PrintRegionTree(string filename) {
            var sw = new System.IO.StreamWriter (filename);
            sw.WriteLine ("Chr,Start,End,Count");
            foreach (var r in tree.PreOrderTraversal()) {
                sw.WriteLine (r.Reference + "," + r.Start + "," + r.End + "," + r.ObservationCount);
            }
            sw.Close ();
        }
        //TODO: Add back in when Bio.Variant namespace is on nuget.
        #if FALSE
        public static Region GetRegionInformation(Variant v) {
            var r = new Region (v.RefName, v.StartPosition, v.EndPosition);
            Region res;
            var exists = tree.TrySearch (r, out res);
            if (!exists)
                return r;
            else
                return res;
        }
        #endif

        #region IDisposable implementation
        public void Dispose ()
        {
            if (bwa != null) {
                bwa.Dispose ();
            }
        }
        #endregion
    }
}

