using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;

using Bio;
using Bio.IO.FastA;
using Bio.Extensions;

namespace Bio.BWA
{
    public class SmallFastaSubSequenceProvider : FastaSubSequenceProvider
    {
        Dictionary<string, byte[]> refData;

        public SmallFastaSubSequenceProvider (string fileName)
        {
            refData = new Dictionary<string, byte[]> ();
            FastAParser fp = new FastAParser ();
            foreach (var seq in fp.Parse (fileName)) {
                var id = seq.ID.Split(' ')[0];
                var seqarr = seq.ToArray ();
                refData [id] = seqarr;
               // var sequence = seq.ConvertToString ();
               // sequence = sequence.ToUpperInvariant ();
               // refData [id] = sequence.Select (z => (byte)z).ToArray ();
            }            
        }
        public override byte[] GetReferenceSection (string refName, long start, long end)
        {
            byte[] data;
            if (!refData.TryGetValue (refName, out data)) {
                throw new ArgumentException ("Reference " + refName + " was not found in the fasta file");
            }
            if (end > data.Length) {
                throw new ArgumentOutOfRangeException ("Requested a location out of range for the reference " + refName);
            }
            if (start > end) {
                throw new ArgumentOutOfRangeException ("start", "Start > End");
            }

            long size = end - start;
            byte[] subsection = new byte[size];
            Array.Copy (data, start, subsection, 0, (int)size);
            return subsection;
        }
    }
}

