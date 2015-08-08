using System;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.IO;
using Bio.IO.FastA;
using Bio;
using System.IO.MemoryMappedFiles;

namespace Bio.IO.FastA
{
    /// <summary>
    /// For large FASTA files, we retrieve subsequences by using a FASTA index and a memory mapped file
    /// instead of grabbing each value directly from the disk.
    /// 
    /// Fasta index.
    /// http://www.htslib.org/doc/samtools.html
    /// </summary>
    public class LargeFastaSubSequenceProvider : FastaSubSequenceProvider
    {
        /// <summary>
        /// A dictionary to look up the index by chromosome position.
        /// </summary>
        Dictionary<string, FastaIndexLine> index;

        /// <summary>
        /// The memory mapped fasta file.
        /// </summary>
        MemoryMappedFile mmf;

        public LargeFastaSubSequenceProvider (string fname, string indexName = null)
        {
            if (indexName == null) {
                indexName = fname + ".fai";
            }
            if (!File.Exists (fname)) {
                throw new FileNotFoundException ("Could not find FASTA file: " + fname);
            } else if(!File.Exists(indexName)) {
                throw new FileNotFoundException ("Could not find fasta index file: " + fname +
                    "\nNote: samtools faidx can be used to generate an index file.");
            }
            index = File.ReadLines (indexName).Select (x => new FastaIndexLine (x)).ToDictionary( z=>z.Name, z=>z);
            try {
                mmf = MemoryMappedFile.CreateFromFile(fname, FileMode.Open);
            }
            catch(Exception thrown) {
                throw new BWA.MEM.BWAException ("Could not open memory mapped file to read fasta.  Error was: \n" +
                thrown.Message + " with stack trace " + thrown.StackTrace +
                "\n\n\nNote that on some operating systems you may need to chmod 777 the file in order to have sufficient permissions to do this.");
            }

        }
        public override byte[] GetReferenceSection(string refName, long start, long end) {
            FastaIndexLine line;
            bool hasValue = index.TryGetValue (refName, out line);
            if (!hasValue) {
                throw new ArgumentException ("Reference: " + refName + " was not found in the index file and cannot be queried");
            }
            if (end > line.Length) {
                throw new ArgumentOutOfRangeException ("Requested a location out of range for the reference " + refName);
            }

            long line_ends_before_start =  start / line.LineLength;
            long line_ends_before_end =  end / line.LineLength;
            long additionalBytesNeeded = line_ends_before_end - line_ends_before_start;

            long size = end - start;
            long file_start = line.Offset + start + start / line.LineLength;
            byte[] seqWithLineEnds = new byte[size + additionalBytesNeeded];
            using (var view = mmf.CreateViewAccessor (file_start, size + additionalBytesNeeded, MemoryMappedFileAccess.Read)) {
                view.ReadArray<byte>(0,seqWithLineEnds, 0, (int)(size + additionalBytesNeeded)); 
            }
            var seq = new byte[size];
            int si = 0;
            for (int i = 0; i < seqWithLineEnds.Length; i++) {
                var bp = seqWithLineEnds [i];
                if (bp != '\n') {
                    // upgrade it to the caps
                    if (bp > 96) {
                        bp = (byte)(bp - 32);
                    }
                    seq [si] = bp;
                    si++;
                }
            }
            return seq;
        }

    }

    /// <summary>
    ///     Fasta index line.
    /// 
    /// 
    /// </summary>
    public class FastaIndexLine {
        public readonly string Name;
        public readonly long Length;
        public readonly long Offset;
        public readonly long LineLength;
        public readonly long BinaryLineLength;

        public FastaIndexLine(string line) {
            var sp = line.Trim ().Split ('\t');
            Name = sp [0];
            Length = Convert.ToInt64 (sp [1]);
            Offset = Convert.ToInt64 (sp [2]);
            LineLength = Convert.ToInt64 (sp [3]);
            BinaryLineLength = Convert.ToInt64 (sp [4]);
            if (BinaryLineLength != LineLength + 1) {
                throw new ArgumentException("This library made assumptions about the binary encoding of the fasta file that have not been met");

            }

        }
    }
}
