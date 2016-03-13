using System;
using System.Runtime.InteropServices;
using System.Linq;
using System.IO;
using Bio.IO.SAM;
using Bio.IO.BAM;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using Bio.Extensions;

namespace Bio.BWA.MEM
{
	public class BWA : IDisposable
	{
		private IntPtr bwaidx; //bwaidx_t*

		//TODO: Convert to the actual structure after profiling to ensure blittable types, right now 
		//the scoring matrix needs to be manually updated to after calls, so this should be wrapped.
		private IntPtr opts;//mem_opt_t*
		
        private bwaidx_t bwaidx_as_struct;

		bntseq_t refSeqs;
		
        private string[] refSeqNames;
		
        /// <summary>
		/// Suffixes of files that should be present if the file is indexed
		/// </summary>
		private string[] requiredFileSuffixes = new string[]{".amb",".ann",".bwt",".pac",".sa"};
		
        /// <summary>
		/// Initializes a new instance of the <see cref="Bio.BWA.MEM.BWA"/> class.  This class will requires
		/// write permissions on the source directory of an index has not yet been made.
		/// </summary>
		/// <param name="fastaFile">Fasta file.</param>
        /// <param name="usePacBioOptions">Use the PacBio options? (equivalent to -x pacbio).</param>
        public BWA (string fastaFile, bool usePacBioOptions=false)
		{
			if (!BitConverter.IsLittleEndian) {
				throw new BWAException("The BWA interface can only be used on little-endian machines");		
				}

			if(String.IsNullOrEmpty(fastaFile) || fastaFile.Contains(" "))
			{
				fastaFile = fastaFile==null ? "File was null" : fastaFile;
				throw new BWAException("The fasta file given to use with BWA either was null or contained spaces, which are not allowed in BWA." +
				                       "\nFasta File was: "+fastaFile);
			}
			loadIndex (fastaFile);
			//load options
			//TODO: allow access to options and changing there-of
            opts = usePacBioOptions ? mem_opt_pacbio_init() : mem_opt_init ();
			if (opts == IntPtr.Zero) {
				throw new BWAException("Could not load the BWA default options");
			}
		}
		
        /// <summary>
		/// Returns the best alignment for that sequence, with information about
		/// other alignments in the meta-data under "NH" meta-tag tag, for number of other alignments.
		/// 
		/// Note that this does not handle paired end data, and like all of bwa mem, the mapping quality does not
		/// take read quality scores in to account.  Also, bwa mem does not places mismatches in cigars, so no assumptions
		/// should be made about what a "60M" means.
		///
		/// </summary>
		/// <returns>The sequence.</returns>
		/// <param name="seq">Sequence to align</param>
		/// <param name="addMetaDataInformat"> Add meta-data? Provides more information at a slight performance cost </param>
		public unsafe SAMAlignedSequence AlignSequence(ISequence seq, bool addMetaDataInformation=true)
		{
			if (seq == null) {
				throw new ArgumentNullException ("seq");
			}
			SAMAlignedSequence toReturn = null;
			//TODO: Validate the alphabet and sequence are of DNA type
			byte[] data = seq.ToArray ();
			fixed(byte * pdata=data)
			{
				//get the alignments
				//TODO: Who is responsible for cleaning out this memory in this ar variable???
				//I think structs do not need to be freed, but best to check
				mem_alnreg_v ar=mem_align1(opts, bwaidx_as_struct.bwt,
				                           bwaidx_as_struct.bns,
				                           bwaidx_as_struct.pac,
				                           data.Length,
				                           (IntPtr)pdata);
				//cycle through them returning the best
				int n_aligns = (int)ar.n;
				//save memory location to free later
				IntPtr initialMemLocation = ar.a;
                int bestScore = int.MinValue;
				for (int i=0; i<n_aligns; i++) {
					mem_alnreg_t curAlign = *(mem_alnreg_t*)ar.a;
                    /* Ignore secondary alignments, but note that the
                     * can be more than one primary alignment.  The criteria is that an alignment is
                     * only secondary if it overlaps with the other alignments query start/end by 50%
                     * (or whatever the mask level option is set to).
                     * 
                     * As a result, I also need to keep track of what the highest scoring PRIMARY alignment is.
                     * In the case of two equivalently scoring primary alignments (yeah, that happens) I select the one with the 
                     * lowest start position.
                     */                      
                    if(curAlign.secondary < 0  && ( 
                        (curAlign.score > bestScore) || (curAlign.score == bestScore && curAlign.rb < toReturn.Pos)))
					{
                        bestScore = curAlign.score;
                        ISequence querySeq = seq;
						//get forward strand positiion and cigar
						//a = mem_reg2aln(opt, idx->bns, idx->pac, ks->seq.l, ks->seq.s, &ar.a[i]); // get forward-strand position and CIGAR
						mem_aln_t a = mem_reg2aln (opts, bwaidx_as_struct.bns,bwaidx_as_struct.pac, data.Length, (IntPtr)pdata, ar.a);
						//unpack the cigar string, which is encoded in the BAM format
						IntPtr originalCigarLocation = a.cigar;
						uint* cigarPointer = (uint*)a.cigar;
						StringBuilder cigarBuilder=new StringBuilder();
						for(int k=0;k<a.n_cigar;++k)
						{
							uint curCigar = *cigarPointer;
							uint opLength = curCigar >> 4;
							char cigarOp = "MIDSH"[(int)(curCigar & 0xf)];
							cigarBuilder.Append (opLength.ToString ());
							cigarBuilder.Append(cigarOp.ToString());
							cigarPointer++;
							 //"0 is +, 1 is -
						}
						//free memory for cigar string
						mem_nd_free_uint (originalCigarLocation);
						//first bit has the reversed or not
                        bool isReversed = (a.isRevAndisAltAndMapQAndNM & 1) == 0 ? false : true;
                        if (isReversed) {
                            querySeq = querySeq.GetReverseComplementedSequence ();
                            querySeq.Metadata.Add (SequenceExtensions.ReversedSequenceMetadataKey, true); 
                        }
						//is_rev:1, is_alt:1, mapq:8, NM:22;
						//next 8 are the mapping quality
                        byte mapq=(byte) ((a.isRevAndisAltAndMapQAndNM >> 2) & 0xFF);
						//and last 23 are edit distance
                        uint editDistance= a.isRevAndisAltAndMapQAndNM >> 10;
						// now let's make a new sequence
						toReturn = new SAMAlignedSequence ();
						// flag includes 0x100 for secondary alignment, and I believe that is the only one it could be at this point
                        toReturn.Flag = (SAMFlags) (a.flag | (isReversed ? (int)SAMFlags.QueryOnReverseStrand : 0)) ;
						toReturn.MapQ = (int)mapq;
						toReturn.CIGAR = cigarBuilder.ToString ();
						toReturn.QName = seq.ID;
                        toReturn.QuerySequence = querySeq;
						toReturn.RName = refSeqNames [a.rid];
						toReturn.Pos = (int)a.pos;
						if (addMetaDataInformation) {
							//Edit distance
							toReturn.Metadata ["NM"] = editDistance;
							//alignment score
							toReturn.Metadata ["AS"] = a.score;
							//other alignments found
							toReturn.Metadata ["NH"] = n_aligns;
							//get second best score here as well
							toReturn.Metadata ["XS"] = a.sub;
						}
					}
					ar.a += sizeof(mem_alnreg_t);
				}
				//free memory
				mem_nd_free_mem_aln_t (initialMemLocation);
			}
			return toReturn;
		}
		private unsafe void loadIndex(string fname)
		{
			//see if we need the index, and load if so
			var filesNotPresent = requiredFileSuffixes.Any (x => !File.Exists (fname + x));
			if (filesNotPresent) {
            	int res=bwa_index (2, new string[] {"index", fname});
				if (res != 0) {
					throw new BWAException ("Could not index fasta file: " + fname);
				}
			}

            //load the index
			bwaidx = _loadIndex (fname);
			if (bwaidx == IntPtr.Zero) {
				throw new BWAException("The BWA files indexed genome failed to load.  Ensure it is in the correct place" +
				                       "and that the bwa index command has been run.");
			}
			bwaidx_as_struct = *(bwaidx_t*)bwaidx;
			getReferenceSequenceInformation ();
		}
		/// <summary>
		/// Gets the reference sequence information by extracting it from the pointers
		/// </summary>
		private unsafe void getReferenceSequenceInformation()
		{
			//the reference sequence is held in an index, the second item of which, bns, is a pointer to a
			// struct of bntseq_t so the trick will be
			//to grab the reference names and lengths from these.
			//mimicing idx->bns->anns[a.rid].name
			refSeqs = *(bntseq_t*)this.bwaidx_as_struct.bns;
			int nRefs = refSeqs.n_seqs;
			if (nRefs == 0) {
				//catastrophic problem here, should never happen.
				this.Dispose (true);
				throw new BWAException ("No reference sequences were loaded by BWA.  Please check your file.");
			}
			refSeqNames = new string[nRefs];
			IntPtr annotations = refSeqs.anns;
			for (int i=0; i<nRefs; i++) {
				bntann1_t refdata = *(bntann1_t*)(annotations);
				//now this is some crap, keep traversing the pointer until we find something null
				List<byte> refName = new List<byte> ();
				sbyte* curSpot = (sbyte*)refdata.name;
				refSeqNames [i] = new string (curSpot);

				annotations += sizeof(bntann1_t);
			}
		}


		/// <summary>
		/// Program taken from c code that updates the scoring matrix in the mem_opts_t structure, this should be called when parameters are changed.
		/// original third argument type was: int8_t mat[25]
		/// </summary>
		/// <param name="a">match score</param>
		/// <param name="b">mimatch penalty</param>
		private void bwa_fill_scmat(int match, int mismatch, sbyte[] mat)
		{
			sbyte a = (sbyte)match;
			sbyte b = (sbyte)mismatch;
			int i, j, k;
			for (i = k = 0; i < 4; ++i) {
				for (j = 0; j < 4; ++j)
					mat[k++] = i == j ? a : (sbyte)-b;
				mat[k++] = (sbyte)-1; // ambiguous base
			}
			for (j = 0; j < 5; ++j) mat[k++] = (sbyte)-1;
		}


		public IntPtr _loadIndex(string filename)
		{
			return bwa_idx_load (filename, 0x7);//load all

		}
		#region CLEANUP_UNMANAGED_RESOURCES
		protected bool disposed=false;
		public void Dispose()
		{
			Dispose (true);
			GC.SuppressFinalize (this);
		}
		protected virtual void Dispose(bool disposing)
		{
			if (!disposed)
			{
				if (disposing) {
					if(opts!=IntPtr.Zero)
					{
					mem_nd_free_opts (opts);
					opts = IntPtr.Zero;
					}
					if(bwaidx!=IntPtr.Zero)
					{
						// Dispose managed resources.
						bwa_idx_destroy (bwaidx);
						bwaidx = IntPtr.Zero;
					}
				}
			}
			disposed = true;
		}
		/// <summary>
		/// Releases unmanaged resources and performs other cleanup operations before the <see cref="Bio.BWA.BWA"/> is
		/// reclaimed by garbage collection.
		/// </summary>
		~BWA()
		{
			this.Dispose (true);
		}
		#endregion

		#region EXTERNAL_FUNCTIONS
		/// <summary>
		/// Function that loads the index file
		/// </summary>
		/// <param name="filenae">Filenae.</param>
		/// <param name="which">Which.</param>
		[DllImport("bwacsharp")]
		private static extern IntPtr bwa_idx_load(string filename, int which); 
		/// <summary>
		/// Create the default options
		/// </summary>
		[DllImport("bwacsharp")]
		private static extern IntPtr mem_opt_init(); 


        /// <summary>
        /// Create the options for PacBio sequencing
        /// This is a later addition to the library, should work but not
        /// well tested yet.
        /// </summary>
        [DllImport("bwacsharp")]
        private static extern IntPtr mem_opt_pacbio_init();
		/// <summary>
	 /// Find the aligned regions for one query sequence
	 ///
	 /// Note that this routine does not generate CIGAR. CIGAR should be
	 /// generated later by mem_reg2aln() below.
	 ///
	 /// @param opt    alignment parameters
	 /// @param bwt    FM-index of the reference sequence
	 /// @param bns    Information of the reference
	 /// @param pac    2-bit encoded reference
	 /// @param l_seq  length of query sequence
	 /// @param seq    query sequence
	 ///
	 /// @return       list of aligned regions.
	/// </summary>
		[DllImport("bwacsharp")]
		internal static extern mem_alnreg_v mem_align1(IntPtr mem_opt, IntPtr bwt, IntPtr bntseq, IntPtr pac, int l_seq, IntPtr seq);
		/**
		 * Generate CIGAR and forward-strand position from alignment region
		 *
		 * @param opt    alignment parameters
		 * @param bns    Information of the reference
		 * @param pac    2-bit encoded reference
		 * @param l_seq  length of query sequence
		 * @param seq    query sequence
		 * @param ar     one alignment region
		 *
		 * @return       CIGAR, strand, mapping quality and forward-strand position
		 */
		[DllImport("bwacsharp")]
		internal static extern mem_aln_t mem_reg2aln(IntPtr mem_opt, IntPtr bns, IntPtr pac, int l_seq, IntPtr seq, IntPtr ptrTo_mem_alnreg_t);
		//cleanup functions
		[DllImport("bwacsharp")]
		internal static extern void mem_nd_free_mem_aln_t(IntPtr pnt);
		[DllImport("bwacsharp")]
		internal static extern void bwa_idx_destroy (IntPtr idx);
		[DllImport("bwacsharp")]
		internal static extern void mem_nd_free_uint(IntPtr cigar);
		[DllImport("bwacsharp")]
		internal static extern void mem_nd_free_opts(IntPtr cigar);
		/// <summary>
		/// The index command, used to create an index file, equivalent to "bwa index filename"
		/// </summary>
		/// <param name="argc">Number of arguments</param>
		/// <param name="argv">string of all the arguments </param>
		[DllImport("bwacsharp")] internal static extern int bwa_index(int argc, string[] argv); // the "index" command
	#endregion
	}
}

