using System;
using System.Linq;
using System.Collections.Generic;
using Bio;
using Bio.Algorithms.Alignment;

namespace ccs_caller
{
    /// <summary>
    /// This class takes alignments
    /// and generates a list of variants.
    /// </summary>
    public static class VariantCaller 
    {
        
        /// <summary>
        /// Given a pairwise sequence alignment, call variants, producing
        /// a list of SNPs and Indels found in the alignment.
        /// </summary>
        /// <param name="alignment"></param>
        /// <returns></returns>
        /// TODO: Passing in the reference separate from the alignment seems bad.... 
        /// .NET Bio used a PairwiseSequenceAlignment class but that seemed too heavy weight, 
        /// alternatively, the  name for these sequences in the shortend sequence meta data
        /// matches the original, so maybe I could use that?  
        public static List<Variant> CallVariants(byte[] refSeq, BPandQV[] querySeq, bool originallyReverseComplemented)
        {
            List<Variant> variants = new List<Variant>();
            // Convert to byte arrays, note we copy as we manipulate them below when
            // left aligning.

            // Validation
            if (refSeq.Length != querySeq.Length) {
                throw new ArgumentException("Alignment passed to variant calling had unequal length sequences");
            }
            validateNoOverlappingGaps(refSeq, querySeq);

            // Left align indels
            leftAlignIndels(refSeq, querySeq);

            // Now call variants.
            var gap = DnaAlphabet.Instance.Gap;
            int i = 0;
            int refPos = 0;
            int queryPos = 0;
            while( i < refSeq.Length)
            {
                if (refSeq[i] == gap)
                {
                    int len = getGapLength(i, refSeq);
                    var bases = getBases(querySeq, i, len);
                    var newVariant = new IndelVariant(refPos - 1, len, bases, IndelType.Insertion, (i == 0 || (i + len) >= refSeq.Length));
                    //newVariant.QV = querySeq [queryPos].QV;
                    variants.Add(newVariant);
                    i += len;
                    queryPos += len;
                }
                else if (querySeq[i].BP == gap)
                {
                    int len = getGapLength(i, querySeq);
                    var bases = getBases(refSeq, i, len);
                    var newVariant = new IndelVariant(refPos - 1, len, bases, IndelType.Deletion, (i == 0 || (i + len) >= refSeq.Length));
                    // TODO: Verify this is the position for this QV value.  It will not be for reverse complement bases
                    //newVariant.QV = querySeq[queryPos+1].QV;
                    variants.Add(newVariant);
                    i += len;
                    refPos += len;
                }
                else
                {
                    if (querySeq[i].BP != refSeq[i])
                    {
                        var newVariant = new SNPVariant(refPos, (char) querySeq[i].BP, (char)refSeq[i]);
                        //newVariant.QV = querySeq[queryPos].QV;
                        variants.Add(newVariant);
                    }
                    i++; refPos++; queryPos++;
                }
            }
            return variants;
        }


        /// <summary>
        /// Given two byte arrays representing a pairwise alignment, shift them so 
        /// that all deletions start as early as possible.  For example:
        /// TTTTAAAATTTT   -> Converts to -> TTTTAAAATTTT
        /// TTTTAA--TTTT                     TTTT--AATTTT
        /// </summary>
        /// <param name="refseq">Reference Sequency</param>
        /// <param name="query">Query Sequence</param>
        /// <returns></returns>
        private static void leftAlignIndels(byte[] refseq, BPandQV[] query)
        {            
            byte gap = DnaAlphabet.Instance.Gap;
            // Keep left aligning until we can't anymore, this is a 
            // do while loop because some downstream left alignments open up
            // further ones upstream, even though this is rare.
            int change_count = 0;
            int loopsThrough = 0;
            do
            {
                loopsThrough++;
                change_count = 0;
                for (int i = 1; i < refseq.Length; i++)
                {
                    if (refseq[i] == gap)
                    {
                        int len = getGapLength(i, refseq);
                        int left_side = i - 1;
                        int right_side = i  - 1 + len;
                        while (left_side >= 0 && refseq[left_side] != gap && (refseq[left_side] == query[right_side].BP))
                        {
                            // Move the gap left.
                            refseq[right_side] = refseq[left_side];
                            refseq[left_side] = gap;
                            left_side--;
                            right_side--;
                            change_count++;
                        }
                        if (loopsThrough > 50) {

                            var s1 = new Sequence(DnaAlphabet.Instance, new string(query.Select(z=>(char)z.BP).ToArray()));
                            var s2 = new Sequence(DnaAlphabet.Instance, refseq);
                            Console.WriteLine(change_count.ToString() +" query");
                            Console.WriteLine(s2.ConvertToString());
                            Console.WriteLine(s1.ConvertToString());
                            throw new Exception();
                        }
                    }
                    else if (query[i].BP == gap)
                    {

                        int len = getGapLength(i, query);
                        int left_side = i - 1;
                        int right_side = i - 1 + len;
                        while (left_side >= 0 && query[left_side].BP != gap && (query[left_side].BP == refseq[right_side]))
                        {
                            // Move the gap left.
                            query[right_side] = query[left_side];
                            query[left_side] = new BPandQV(gap, 0);
                            left_side--;
                            right_side--;
                            change_count++;
                        }
                        if (loopsThrough > 50) {

                            var s1 = new Sequence(DnaAlphabet.Instance, new string(query.Select(z=>(char)z.BP).ToArray()));
                            var s2 = new Sequence(DnaAlphabet.Instance, refseq);
                            Console.WriteLine(change_count.ToString() +" query");
                            Console.WriteLine(s2.ConvertToString());
                            Console.WriteLine(s1.ConvertToString());
                            throw new Exception();
                        }
                    }
                }
            } while (change_count > 0);
        }

        /// <summary>
        /// Simple check that the alignment does not have a gap on top of a
        /// gap, which violates several assumptions.
        /// </summary>
        /// <param name="seq1"></param>
        /// <param name="seq2"></param>
        private static void validateNoOverlappingGaps(byte[] seq1, BPandQV[] seq2)
        {
            var gap = DnaAlphabet.Instance.Gap;
            for(int i=0;i<seq1.Length;i++)
            {
                if (seq1[i] == gap && seq2[i].BP == gap)
                    throw new Exception("You have an alignment with overlapping gaps.  Input problem!");
            }
        }

        /// <summary>
        /// Given the start position of a gap, returns how long it is.
        /// For example:
        /// 
        /// AAAA---TTTT returns 3.
        /// </summary>
        /// <param name="pos">0 indexed</param>
        /// <param name="array"></param>
        /// <returns></returns>
        private static int getGapLength(int pos, byte[] array)
        {
            var gap = DnaAlphabet.Instance.Gap;
            int len =1;
            while (pos <= array.Length)
            {
                if (array[++pos]==gap)
                {
                    len += 1;
                }
                else
                {
                    break;
                }
            }
            return len;
        }

        private static int getGapLength(int pos, BPandQV[] array)
        {
            var gap = DnaAlphabet.Instance.Gap;
            int len =1;
            while (pos <= array.Length)
            {
                if (array[++pos].BP == gap)
                {
                    len += 1;
                }
                else
                {
                    break;
                }
            }
            return len;
        }

        /// <summary>
        /// Converts a subset of bases in the array into a string.
        /// </summary>
        /// <param name="array"></param>
        /// <param name="position"></param>
        /// <param name="length"></param>
        /// <returns></returns>
        private static string getBases(byte[] array, int position, int length)
        {
            char[] chars = new char[length];
            for(int i=0; i<length; i++)
            {
                chars[i] = (char)array[i+position];
            }
            return new string(chars);
        }
        private static string getBases(BPandQV[] array, int position, int length)
        {
            char[] chars = new char[length];
            for(int i=0; i<length; i++)
            {
                chars[i] = (char)array[i+position].BP;
            }
            return new string(chars);
        }
    }
}

