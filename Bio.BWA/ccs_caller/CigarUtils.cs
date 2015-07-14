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

