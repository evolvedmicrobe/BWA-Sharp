using System;

namespace Bio.IO.SAM
{
    public static class CigarOperations
    {
        /// <summary>
        /// Alignment position match, but sequence can be match or mismatch
        /// </summary>
        public const char ALN_MATCH = 'M';

        /// <summary>
        /// Insertion to the reference
        /// </summary>
        public const char INSERTION = 'I';

        /// <summary>
        /// Deletion from the reference
        /// </summary>
        public const char DELETION = 'D';

        /// <summary>
        /// Skipped region from the reference
        /// </summary>
        public const char SKIPPED = 'N';

        /// <summary>
        /// Soft Clipping, clipped sequences present in sequence
        /// </summary>
        public const char SOFT_CLIP = 'S';

        /// <summary>
        /// Hard clipping, clipped sequences not present in SEQ.
        /// </summary>
        public const char HARD_CLIP = 'H';

        /// <summary>
        /// Padding (silent deletion from padded reference)
        /// </summary>
        public const char PADDING = 'P';

        /// <summary>
        /// Sequence match
        /// </summary>
        public const char SEQ_MATCH = '=';

        /// <summary>
        /// Sequence mismatch
        /// </summary>
        public const char SEQ_MISMATCH = 'X';


        public static bool OperationIsMatch(char operation)
        {
            return operation == SEQ_MATCH || operation == SEQ_MISMATCH || operation == ALN_MATCH;
        }
    }
}
