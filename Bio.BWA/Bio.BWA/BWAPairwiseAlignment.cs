using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;

using Bio.Algorithms.Alignment;
using Bio.BWA;
using Bio.IO;
using Bio.IO.SAM;
using Bio.Extensions;
namespace Bio.BWA
{
    public class BWAPairwiseAlignment : IPairwiseSequenceAlignment
    {
        public string Reference {
            get {
                return AlignedSAMSequence.RName;
            }
        }

        /// <summary>
        /// 0 based start of alignment on reference.
        /// </summary>
        /// <value>The start.</value>
        public long Start {
            get{ return AlignedSAMSequence.Pos; }
        }
        /// <summary>
        /// 0 based exclusive end of the alignment (+1 last aligned position).
        /// </summary>
        /// <value>The end.</value>
        public long End {
            get { return this.AlignedSAMSequence.RefEndPos + 1; }
        }
        public readonly bool OriginallyReverseComplemented;

        public int MapQV {
            get { return this.AlignedSAMSequence.MapQ; }
        }

        public Sequence AlignedRefSeq;
        public ISequence AlignedQuerySeq;
        PairwiseAlignedSequence pas;
        public SAMAlignedSequence AlignedSAMSequence;

        public BWAPairwiseAlignment (SAMAlignedSequence aln, Sequence alignedRefSeq, ISequence alignedQualitySequence)
        {
            if (aln == null) {
                throw new ArgumentNullException ("aln");
            } else if (alignedRefSeq == null) {
                throw new ArgumentNullException ("alignedRefSeq");
            } else if (alignedQualitySequence == null) {
                throw new ArgumentNullException ("alignedQualitySequence");
            }
            AlignedSAMSequence = aln;
            this.OriginallyReverseComplemented = ((aln.Flag & SAMFlags.QueryOnReverseStrand) == SAMFlags.QueryOnReverseStrand);
            if (OriginallyReverseComplemented) {
                alignedQualitySequence.MarkAsReverseComplement ();
            }
            this.AlignedRefSeq = alignedRefSeq;
            this.AlignedQuerySeq = alignedQualitySequence;

            pas = new PairwiseAlignedSequence ();
            pas.FirstSequence = alignedRefSeq;
            pas.SecondSequence = alignedQualitySequence;

        }

        #region IPairwiseSequenceAlignment implementation

        public void AddSequence (PairwiseAlignedSequence pairwiseAlignedSequence)
        {
            throw new NotImplementedException ();
        }

        public System.Collections.Generic.IList<PairwiseAlignedSequence> PairwiseAlignedSequences {
            get {
                return new List<PairwiseAlignedSequence> () { pas };
            }
        }


        public ISequence FirstSequence {
            get {
                throw new NotImplementedException ();
                }
        }

        public ISequence SecondSequence {
            get {
                return this.AlignedSAMSequence.QuerySequence;
            }
        }

        #endregion

        #region ICollection implementation

        public void Add (PairwiseAlignedSequence item)
        {
            throw new NotImplementedException ();
        }

        public void Clear ()
        {
            throw new NotImplementedException ();
        }

        public bool Contains (PairwiseAlignedSequence item)
        {
            throw new NotImplementedException ();
        }

        public void CopyTo (PairwiseAlignedSequence[] array, int arrayIndex)
        {
            throw new NotImplementedException ();
        }

        public bool Remove (PairwiseAlignedSequence item)
        {
            throw new NotImplementedException ();
        }

        public int Count {
            get {
                return 2;
            }
        }

        public bool IsReadOnly {
            get {
                throw new NotImplementedException ();
            }
        }

        #endregion

        #region IEnumerable implementation

        public System.Collections.Generic.IEnumerator<PairwiseAlignedSequence> GetEnumerator ()
        {
            throw new NotImplementedException ();
        }

        #endregion

        #region IEnumerable implementation

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator ()
        {
            throw new NotImplementedException ();
        }

        #endregion

        #region ISequenceAlignment implementation

        public System.Collections.Generic.IList<IAlignedSequence> AlignedSequences {
            get {
                throw new NotImplementedException ();
            }
        }

        public System.Collections.Generic.IList<ISequence> Sequences {
            get {
                throw new NotImplementedException ();
            }
        }

        public System.Collections.Generic.IDictionary<string, object> Metadata {
            get {
                throw new NotImplementedException ();
            }
        }

        public object Documentation {
            get {
                throw new NotImplementedException ();
            }
            set {
                throw new NotImplementedException ();
            }
        }

        #endregion
    }
}

