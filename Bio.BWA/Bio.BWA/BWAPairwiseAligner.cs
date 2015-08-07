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

        public BWAPairwiseAligner (string fasta_name, bool recordCoverageIntervals=false)
        {
            this.fasta_name = fasta_name;
            this.recordCoverageIntervals = recordCoverageIntervals;
            if (recordCoverageIntervals) {
                tree = new RegionTree ();
            }
            bwa = new Bio.BWA.MEM.BWA (fasta_name);

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
            elements = CigarUtils.FilterLeadingAndTrailingInsertions (elements);

            int length = CigarUtils.GetAlignmentLengthAfterClipping (elements);
            var q = new List<BPandQV>(length);
            var r = new List<byte>(length);

            int curRef = 0;
            int curQuery = 0; //location on query
            var qseq = aln.QuerySequence as IQualitativeSequence;
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
                        var qv = isQual ? qseq.GetEncodedQualityScore (curQuery) : (byte)0;
                        q.Add (new BPandQV (bp, qv));
                        curQuery++;
                    } else {
                        q.Add (new BPandQV ((byte)'-', 0));
                    }
                }

            }
            return new Tuple<List<byte>, List<BPandQV>> (r, q);
        }

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

