using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace Bio.IO.SAM
{
    public struct CigarElement
    {
        public int Length;
        public char Operation;

        public CigarElement (char operation, int length)
        {
            this.Length = length;
            this.Operation = operation;
        }
    }

    public static class CigarUtils
    {
        
        public static List<CigarElement> GetCigarElements(string cigar)
        {
            var elements = new List<CigarElement> (7);
            int start = 0;
            for (int i = 1; i < cigar.Length; i++) {
                char ch = cigar [i];
                if (!Char.IsDigit (ch)) {
                    var str_len = i - start;
                    var length = int.Parse (cigar.Substring (start, str_len), CultureInfo.InvariantCulture);
                    elements.Add (new CigarElement (ch, length));
                    start = ++i;//increment i as the next element will be a number
                }
            }
            return elements;
        }

        /// <summary>
        /// Apparently BWA is quite content to return a CIGAR like "6S1I20M1D50M" 
        /// where there is both soft clipping and an insertion at the start. To avoid this, I will filter out any
        /// leading insertions.  This appears to be a rare event.
        /// 
        /// Not sure if this can also happen with deletions.
        /// </summary>
        /// <returns>The trailing insertions.</returns>
        public static List<CigarElement> FilterLeadingAndTrailingInsertions(List<CigarElement> cigars) {
            
            // Filter leading
            while (cigars [0].Operation == CigarOperations.SOFT_CLIP && 
                cigars.Count > 1 &&
                cigars [1].Operation == CigarOperations.INSERTION) 
            {
                var old = cigars [0];
                old.Length += cigars [1].Length;
                cigars [0] = old;
                cigars.RemoveAt (1);
            }

            // Filter trailing
            while (cigars [cigars.Count - 1].Operation == CigarOperations.SOFT_CLIP && 
                cigars.Count > 1 &&
                cigars [cigars.Count -2].Operation == CigarOperations.INSERTION) 
            {
                var old = cigars [cigars.Count - 1];
                old.Length += cigars [cigars.Count - 2].Length;
                cigars [cigars.Count - 1] = old;
                cigars.RemoveAt (cigars.Count - 2);
            }

            if (cigars [0].Operation == CigarOperations.SOFT_CLIP &&
                cigars.Count > 1 &&
                cigars [1].Operation == CigarOperations.DELETION) {
                throw new FormatException ("BWA return a cigar string with a softclip followed by a deletion.");
            }
            if (cigars [cigars.Count - 1].Operation == CigarOperations.SOFT_CLIP &&
                cigars.Count > 1 &&
                cigars [cigars.Count - 2].Operation == CigarOperations.DELETION) {
                throw new FormatException ("BWA return a cigar string with a deletion followed by a softclip.");
            }
            return cigars;
        }
        public static int GetAlignmentLengthAfterClipping(IEnumerable<CigarElement> elements)
        {            
            int length = 0;
            foreach (var e in elements) {
                char type = e.Operation;
                if (type == CigarOperations.PADDING ||
                    type == CigarOperations.HARD_CLIP ||
                    type == CigarOperations.SKIPPED) {
                    throw new FormatException ("Padding, hard clipping and skipped bases are not currently supported.");
                } else if (type != CigarOperations.SOFT_CLIP) {
                    length += e.Length;
                }
            }
            return length;
        }

        public static string CreateCigarString(IEnumerable<CigarElement> elements) {
            return String.Join("",elements.Select(x=> x.Length.ToString()+x.Operation.ToString()));
        }
        /// <summary>
        /// These are the CIGAR elements "MDNX=" can change the alignment length on the reference
        /// this check if it is one of these elements.
        /// </summary>
        /// <returns><c>true</c>, if elementis_ MDN x_ equal was cigared, <c>false</c> otherwise.</returns>
        /// <param name="element">Element.</param>
        public static bool CigarElementis_MDNX_Equal(char element)
        {
            //"MDNX=";
            return element == CigarOperations.ALN_MATCH || element == CigarOperations.DELETION || element == CigarOperations.SKIPPED ||
                element == CigarOperations.SEQ_MISMATCH || element == CigarOperations.SEQ_MATCH;
        }

        public static bool NoInformationCigar(string cigar) {
            return String.IsNullOrEmpty (cigar) || cigar == "*";
        }
    }
}
