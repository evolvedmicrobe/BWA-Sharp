using System;
using System.Runtime.InteropServices;
namespace Bio.BWA
{


	//BWA_IDX_ALL= 0x7

	//original below
//	typedef struct {
//		bwtint_t primary; // S^{-1}(0), or the primary index of BWT
//		bwtint_t L2[5]; // C(), cumulative count
//		bwtint_t seq_len; // sequence length
//		bwtint_t bwt_size; // size of bwt, about seq_len/4
//		uint32_t *bwt; // BWT
//		// occurance array, separated to two parts
//		uint32_t cnt_table[256];
//		// suffix array
//		int sa_intv;
//		bwtint_t n_sa;
//		bwtint_t *sa;
//	} bwt_t;

	using bwtint_t= System.UInt64;//originall unsigned char
	[StructLayout(LayoutKind.Sequential)]
	internal struct bwt_t
	{
		bwtint_t primary; // S^{-1}(0), or the primary index of BWT
		//	bwtint_t L2[5]; // C(), cumulative count
		internal bwtint_t L2_0;
		internal bwtint_t L2_1;
		internal bwtint_t L2_2;
		internal bwtint_t L2_3;
		internal bwtint_t L2_4;
		internal bwtint_t seq_len; // sequence length
		internal bwtint_t bwt_size; // size of bwt, about seq_len/4
		//this seems wrong to me, why uint32 for a pointer????
		uint bwt_pointer;
		//uint32_t *bwt; // BWT
		// occurance array, separated to two parts
		//uint32_t cnt_table[256];
		[MarshalAs(UnmanagedType.ByValArray, SizeConst=256)] 
		public uint[] cnt_table;
		// suffix array
		int sa_intv;
		//bwtint_t *sa;
		bwtint_t n_sa;
	}


//	typedef struct {
//		bwt_t    *bwt; // FM-index
//		bntseq_t *bns; // information on the reference sequences
//		uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base
//	} bwaidx_t;
	[StructLayout(LayoutKind.Sequential)]
	internal struct bwaidx_t {
			//bwt_t    *bwt; // FM-index
		/// <summary>
		///  FM-index
		/// </summary>
		internal IntPtr bwt;
		/// <summary>
		///  information on the reference sequences
		/// </summary>
		internal IntPtr bns;
		/// <summary>
		/// the actual 2-bit encoded reference sequences with 'N' converted to a random base
		/// </summary>
		internal IntPtr pac;
	} 


//	typedef struct {
//		int64_t rb, re; // [rb,re): reference sequence in the alignment
//		int qb, qe;     // [qb,qe): query sequence in the alignment
//		int score;      // best local SW score
//		int truesc;     // actual score corresponding to the aligned region; possibly smaller than $score
//		int sub;        // 2nd best SW score
//		int csub;       // SW score of a tandem hit
//		int sub_n;      // approximate number of suboptimal hits
//		int w;          // actual band width used in extension
//		int seedcov;    // length of regions coverged by seeds
//		int secondary;  // index of the parent hit shadowing the current hit; <0 if primary
                  
	[StructLayout(LayoutKind.Sequential)]
	internal struct mem_alnreg_t {
		/// <summary>
		/// Reference begin
		/// </summary>
		internal long rb;
		/// <summary>
		/// Reference end
		/// </summary>
		internal long re; // [rb,re): reference sequence in the alignment
			/// <summary>
			/// Query begin
			/// </summary>
		internal int qb;
		/// <summary>
		/// Query end
		/// </summary>
		internal int qe;     // [qb,qe): query sequence in the alignment
			/// <summary>
			/// best local SW score
			/// </summary>
		internal int score;      // best local SW score
		/// <summary>
		/// actual score corresponding to the aligned region
		/// </summary>
		internal int truesc;     // actual score corresponding to the aligned region; possibly smaller than $score
	/// <summary>
	/// 2nd best SW score
	/// </summary>
		internal int sub;        // 2nd best SW score
		/// <summary>
		/// SW score of a tandem hit/
		/// </summary>
		internal int csub;       
		/// <summary>
		/// Approximate number of suboptimal hits
		/// </summary>
		internal int sub_n;      // approximate number of suboptimal hits
		/// <summary>
		/// Actual band width used in extension
		/// </summary>
		internal int w;          // actual band width used in extension
		/// <summary>
		/// Length of regions coverged by seeds
		/// </summary>
		internal int seedcov;    // length of regions coverged by seeds
		/// <summary>
		/// ndex of the parent hit shadowing the current hit; <0 if primary
		/// </summary>
		internal int secondary;  // index of the parent hit shadowing the current hit; <0 if primary
   	}
	//typedef struct { size_t n, m; mem_alnreg_t *a; } mem_alnreg_v;
	[StructLayout(LayoutKind.Sequential)]
	internal struct mem_alnreg_v
	{
		internal UIntPtr n,m;
		/// <summary>
		/// Points to mem_aln_reg_type
		/// </summary>
		internal IntPtr a;

	}

//	typedef struct { // This struct is only used for the convenience of API.
//		int64_t pos;     // forward strand 5'-end mapping position
//		int rid;         // reference sequence index in bntseq_t; <0 for unmapped
//		int flag;        // 
//		uint32_t is_rev:1, mapq:8, NM:23; // is_rev: whether on the reverse strand; mapq: mapping quality; NM: edit distance
//		int n_cigar;     // number of CIGAR operations
//		uint32_t *cigar; // CIGAR in the BAM encoding: opLen<<4|op; op to integer mapping: MIDSH=>01234
//
//		int score, sub;
//	} mem_aln_t;
	[StructLayout(LayoutKind.Sequential)]
	internal struct mem_aln_t
	{
		/// <summary>
		/// forward strand 5'-end mapping position
		/// </summary>
		internal long pos;
		/// <summary>
		/// reference sequence index in bntseq_t; <0 for unmapped
		/// </summary>
		internal int rid;
		/// <summary>
		/// Extra flag
		/// </summary>
		internal int flag;
		/// <summary>
		/// is_rev:1, mapq:8, NM:23;
		/// </summary>
		internal uint isRevAndMapQAndNM;
		/// <summary>
		/// number of cigar operations
		/// </summary>
		internal int n_cigar;
		/// <summary>
		/// Pointer to cigar in bam encoding
		/// </summary>
		internal IntPtr cigar;

		internal int score;
		internal int sub;
	}


	//	#define __KSEQ_TYPE(type_t)						\
//typedef struct {							\
//kstring_t name, comment, seq, qual;		\
//int last_char;							\
//kstream_t *f;							\
//} kseq_t;



}

