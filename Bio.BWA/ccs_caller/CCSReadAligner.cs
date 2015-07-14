using System;
using Bio;
using Bio.IO.FastQ;
using Bio.BWA;
using Bio.BWA.MEM;
using Bio.IO.SAM;
using Bio.IO.FastA;

namespace ccs_caller
{
    
    public static class CCSReadAligner
    {
        static Bio.BWA.MEM.BWA bwa = new BWA("/Users/nigel/BroadBundle/human_g1k_v37.fasta"); 
        const string fname = "/Users/nigel/BroadBundle/human_g1k_v37.fasta";
        static RegionTree tree = new RegionTree();
        internal static FastaIndex fi = new FastaIndex (fname);
        static int UnAlignedReads =0;
        static int AlignedReads = 0;
        public static AlignedCCSRead AlignRead(IQualitativeSequence read, int numPasses)
        {
            var aln = bwa.AlignSequence (read);
            if (aln == null) {
                UnAlignedReads++;
                return null;
            }
            else {
                AlignedReads++;

                var refSeq = fi.GetReferenceSection (aln.RName, aln.Pos, aln.RefEndPos + 1);
                var res = new AlignedCCSRead (aln, refSeq, numPasses);
                var reg = new Region (aln.RName, aln.Pos, aln.RefEndPos);
                tree.Add (reg);
                return res;
            }
        }
        public static void PrintRegionTree(string filename) {
            var sw = new System.IO.StreamWriter (filename);
            sw.WriteLine ("Chr,Start,End,Count");
            foreach (var r in tree.PreOrderTraversal()) {
                sw.WriteLine (r.Reference + "," + r.Start + "," + r.End + "," + r.ObservationCount);
            }
            sw.Close ();
        }
        public static Region GetRegionInformation(Variant v) {
            var r = new Region (v.RefName, v.StartPosition, v.EndPosition);
            Region res;
            var exists = tree.TrySearch (r, out res);
            if (!exists)
                return r;
            else
                return res;
        }
    }
}

