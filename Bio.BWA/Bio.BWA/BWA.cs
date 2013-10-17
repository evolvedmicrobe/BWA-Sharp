using System;
using System.Runtime.InteropServices;

namespace Bio.BWA
{
	public class BWA
	{
		private IntPtr bwaidx; //bwaidx_t*
		private IntPtr opt;//mem_opt_t*
		public unsafe BWA ()
		{
			if (!BitConverter.IsLittleEndian) {
				throw new Exception("We assumed little endian!!!");		
				}
			Console.WriteLine (System.IO.Directory.GetCurrentDirectory ());
			bwaidx = LoadIndex ("/Users/ndelaney/bwa-0.7.5a/TestData/MT.fasta");
			if (bwaidx == IntPtr.Zero) {
				throw new BWAException("The BWA files indexed genome failed to load.  Ensure it is in the correct place" +
					"and that the bwa index command has been run.");
			}
			opt = mem_opt_init ();
			bwaidx_t tmp = *(bwaidx_t*)bwaidx;
			///from position 5110 and of length 65, the first 60 bases match, the next 5 do not.
			string test = "CGCATTCCTACTACTCAACTTAAACTCCAGCACCACGACCCTACTACTATCTCGCACCTGCTTTT";
			//it's output should be equal to:"test	+	MT	5110	60	60M5S	0"
			byte[] data = System.Text.Encoding.ASCII.GetBytes (test);
			fixed(byte * pdata=data)
			{
				mem_alnreg_v ar=mem_align1(opt,tmp.bwt,tmp.bns,tmp.pac,data.Length,(IntPtr)pdata);
			
				int n_aligns = (int)ar.n;
				for (int i=0; i<n_aligns; i++) {
					mem_alnreg_t curAlign = *(mem_alnreg_t*)ar.a;
					if(curAlign.secondary<0)//ignore secondary alignments
					{
						//a = mem_reg2aln(opt, idx->bns, idx->pac, ks->seq.l, ks->seq.s, &ar.a[i]); // get forward-strand position and CIGAR
						mem_aln_t a = mem_reg2aln (opt, tmp.bns, tmp.pac, data.Length, (IntPtr)pdata, ar.a);

					}
					ar.a += sizeof(mem_alnreg_t);
					Console.WriteLine (curAlign.rb);
				}
			}
				Console.WriteLine ("Success!");
		}
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

		public IntPtr LoadIndex(string filename)
		{
			return bwa_idx_load (filename, 0x7);//load all

		}
	

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
	}

}

