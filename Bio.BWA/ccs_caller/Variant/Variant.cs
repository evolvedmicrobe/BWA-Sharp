using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio;

namespace ccs_caller
{
    /// <summary>
    /// Top level class for holding variant information.  
    /// This is implemented in two sub classes, SNPVariant and IndelVariant.
    /// </summary>
    public abstract class Variant : IComparable<Variant>, IEquatable<Variant>
    {

        public string RefName;

        /// <summary>
        /// SNP, indel, etc.
        /// Now redundant with the Variant Type
        /// </summary>
        public VariantType Type { get; protected set; }

        /// <summary>
        /// 0-based start position of variant.  For Indels, this is the left most position 
        /// BEFORE the event. e.g.
        /// 
        /// AAATTTAAA   -> is a deletion of length 3 starting at position 2.
        /// AAA---AAA
        /// 
        /// A--TTAAA -> is an insertion of length 2 starting at position 0.
        /// ACCTTAAA
        /// </summary>
        public int StartPosition { get; set; }

        /// <summary>
        /// O-based end index (same as start for SNPs).
        /// A SNP is of length 1, a deletion of length 2 is 2, etc.
        /// </summary>
        public int Length { get; protected set; }

        /// <summary>
        /// The position in the alignment where this variant ends.
        /// For a SNP, this is the position AFTER the SNP.  Same for indels.
        /// </summary>
        public int EndPosition
        {
            get { return StartPosition + Length; }
        }

        /// <summary>
        /// Is the variant at the very start or end of an alignment? 
        /// (that is it was called based on the first or last base seen on 
        /// either sequence in the alignment.)
        /// These can have certain pathologies so we take note and keep an eye on them.
        /// They should almost always be excluded by the alignment algorithm clipping at the ends.
        /// </summary>
        public bool AtEndOfAlignment { get; protected set; }

        /// <summary>
        /// The QV value for this call, if it exists.
        /// </summary>
        /// <value>The Q.</value>
        public int QV {get; set;}
       

        public Variant(int position,  bool atAlignmentEnd = false)
        {
            StartPosition = position;
            AtEndOfAlignment = atAlignmentEnd;
        }

        public override int GetHashCode ()
        {
            var h1 = RefName.GetHashCode ();
            var h2 = StartPosition.GetHashCode ();
            var h3 = EndPosition.GetHashCode ();
            return h1 ^ h2 ^ h3;
        }

        public override bool Equals (object obj)
        {
            var res = obj as Variant;
            if (res != null)
                return this.Equals (res);
            return false;
        }
        #region IComparable implementation

        public int CompareTo (Variant other)
        {
            var nameComp = String.CompareOrdinal(RefName,other.RefName);
            if (nameComp == 0) {
                return StartPosition.CompareTo (other.StartPosition);
            }
            return nameComp;
        }

        #endregion

        #region IEquatable implementation

        bool IEquatable<Variant>.Equals (Variant other)
        {
            var otherComp = CompareTo (other);
            if (otherComp == 0) {
                return this.StartPosition == other.StartPosition && this.Type == other.Type && this.EndPosition == other.EndPosition;
            }
            return false;
        }

        #endregion
    }

    public enum VariantType
    {
        SNP,
        INDEL,
        Complex
    }
}
