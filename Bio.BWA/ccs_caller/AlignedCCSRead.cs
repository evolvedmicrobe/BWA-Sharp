using System;
using Bio;
using Bio.IO.FastQ;
using Bio.BWA;
using Bio.BWA.MEM;
using Bio.IO.SAM;
using Bio.IO.FastA;
using System.Linq;
using System.Collections.Generic;
using Bio.Extensions;
using System.Diagnostics;


namespace ccs_caller
{
    [DebuggerDisplay("{Reference} : {Start}-{End}")]
    public class AlignedCCSRead : IComparable<AlignedCCSRead>
    {
        public readonly string Reference;
        public readonly long Start;
        public readonly long End;
        public readonly bool OriginallyReverseComplemented;
        public byte[] AlignedReferenceSeq;
        public BPandQV[] AlignedQuerySeq;
        public int MapQV;
        public ZMWInfo Zmw;

        public readonly List<Variant> Variants;

        //readonly SAMAlignedSequence aln;


        public AlignedCCSRead (SAMAlignedSequence aln, byte[] refSeq, int numPasses)
        {
            //this.aln = aln;
            this.Start = aln.Pos;
            this.End = aln.RefEndPos + 1;
            this.Reference = aln.RName;
            this.MapQV = aln.MapQ;
            this.OriginallyReverseComplemented = ((aln.Flag & SAMFlags.QueryOnReverseStrand) == SAMFlags.QueryOnReverseStrand);


            expandAlignment (refSeq, aln);
            this.Variants = VariantCaller.CallVariants (AlignedReferenceSeq, AlignedQuerySeq, OriginallyReverseComplemented);
            foreach (var v in Variants) {
                v.StartPosition = v.StartPosition + aln.Pos;
                v.RefName = aln.RName;

                var q = v as IndelVariant;
                if (q != null && numPasses>50) {
//                    Console.WriteLine ("\nNewVariant");
//                    Console.WriteLine (q.StartPosition);
//                    Console.WriteLine (numPasses);
//                    Console.WriteLine (q.InHomopolymer);
//                    Console.WriteLine (q.HomopolymerLengthInReference);
//                    Console.WriteLine (q.HomopolymerBase + " " + q.QV);
//                    Console.WriteLine (new String (CCSReadAligner.fi.GetReferenceSection (q.RefName, q.StartPosition, q.StartPosition + 50).Select (x => (char)x).ToArray ()));
//                    Console.WriteLine ("Q: " + new String( AlignedQuerySeq.Select (z => (char)z.BP).ToArray ()));
//                     Console.WriteLine ("R: "+ new String( AlignedReferenceSeq.Select (z => (char)z).ToArray ()));
//                    Console.WriteLine (q.StartPosition);
                }
            }

        }

        void expandAlignment(byte[] refSeq, SAMAlignedSequence aln)
        {
            //Console.WriteLine(new String(refSeq.Select(z=>(char)z).ToArray()));
           // Console.WriteLine(new String(this.aln.QuerySequence.Select(z=>(char)z).ToArray()));
                
            List<BPandQV> q = new List<BPandQV> (refSeq.Length + 10);
            List<byte> r = new List<byte> (refSeq.Length + 10);

            if (string.IsNullOrWhiteSpace (aln.CIGAR) || aln.CIGAR.Equals ("*")) {
                throw new Exception ();
            }
            int curRef = 0;
            int curQuery = 0;//location on query
            var qseq = aln.QuerySequence as QualitativeSequence;
            var qstr = qseq.ConvertToString ();
            var elements = CigarUtils.GetCigarElements (aln.CIGAR);

            // Account for soft clipping in offset
            foreach (var e in elements) {
                char ch = e.Operation;
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
                        var qv = qseq.GetEncodedQualityScore (curQuery);
                        q.Add (new BPandQV (bp, qv));
                        curQuery++;
                    } else {
                        q.Add (new BPandQV ((byte)'-', 0));
                    }
                }

            }
            this.AlignedQuerySeq = q.ToArray ();
            this.AlignedReferenceSeq = r.ToArray ();
            //Console.WriteLine ("Q: " + new String( AlignedQuerySeq.Select (z => (char)z.BP).ToArray ()));
           /// Console.WriteLine ("R: "+ new String( AlignedReferenceSeq.Select (z => (char)z).ToArray ()));

        }
        #if FALSE
        public int? GetReadIndexOfAlignmentPosition (int referencePosition)
        {
            if (string.IsNullOrWhiteSpace (aln.CIGAR) || aln.CIGAR.Equals ("*")) {
                return null;
            }
            long curRef = Start;//location on ref
            int curQuery = 0;//location on query
            if (curRef == referencePosition) {
                return curQuery;
            }
            var elements = CigarUtils.GetCigarElements (aln.CIGAR);
            string CIGARforClen = "MDNX=";
            foreach (var e in elements) {
                char ch = e.Operation;
                if (ch == 'P' || ch == 'N') {
                    throw new Exception ("Not built to handle clipping yet");
                }
                int start = 0;
                int end = 0;
                int len;
                if (CigarUtils.CigarElementis_MDNX_Equal(ch)) {
                    len = e.Length;
                    bool addRef = false;
                    bool addQuery = false;
                    if (ch == 'M' || ch == '=' || ch == 'X') {
                        addRef = true;
                        addQuery = true;
                    } else if (ch == 'D') {
                        addRef = true;
                    } else if (ch == 'I' || ch == 'S' || ch == 'H') {
                        addQuery = true;
                    }
                    for (int j = 0; j < len; j++) {
                        if (addRef)
                            curRef++;
                        if (addQuery)
                            curQuery++;
                        if (curRef == referencePosition && addQuery) {
                            return curQuery;
                        }
                    }
                    if (curRef > referencePosition) {
                        return null;
                    }
                }
            }
            return null;
        }
        #endif

        #region IComparable implementation

        public int CompareTo (AlignedCCSRead other)
        {
            var cmp1 = String.CompareOrdinal(Reference, other.Reference);
            if (cmp1 == 0) {
                return End.CompareTo (other.End);
            }
            return cmp1;
        }

        #endregion


    }
}

