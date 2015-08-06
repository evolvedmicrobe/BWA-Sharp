using System;

namespace Bio.IO.FastA
{
    public abstract class FastaSubSequenceProvider
    {
        abstract public byte[] GetReferenceSection(string refName, long start, long end); 
    }
}

