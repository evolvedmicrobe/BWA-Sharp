using System;

namespace ccs_caller
{
    public class Region : IComparable<Region>
    {
        public string Reference;
        public long Start;
        public long End;
        public long ObservationCount;

        public Region (string reference, long start, long end)
        {
            if (start > end) {
                throw new ArgumentException ("start was after end");
            }
            Reference = reference;
            Start = start;
            End = end;
            ObservationCount = 1;
        }
        public void Merge(Region other) {
            Start = Math.Min(Start, other.Start);
            End = Math.Max (End, other.End);
            ObservationCount++;
        }
        public override bool Equals (object obj)
        {
            var q = obj as Region;
            if (q != null) {
                return (this.CompareTo (q) == 0);
            }
            return base.Equals (obj);
        }
        #region IComparable implementation
        /// <Docs>Compare two intervals</Docs>
        /// <para>Returns 0 if the intervals overlap, negative if less, postive if more.</para>
        /// <summary>
        /// Compares to.
        /// </summary>
        /// <returns>The to.</returns>
        /// <param name="other">Other.</param>
        public int CompareTo (Region other)
        {
            var res = Reference.CompareTo (other.Reference);
            if (res == 0) {
                if (End > other.Start && Start < other.End) {
                    return 0;
                }
                return Start.CompareTo (other.Start);
            }
            return res;
        }
        #endregion
    }
}

