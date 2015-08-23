using System;
using System.Runtime.InteropServices;

namespace Bio.BWA.MEM
{
	//from: typedef uint64_t bwtint_t;
	using bwtint_t= System.UInt64;//originall unsigned char


	//BWA_IDX_ALL= 0x7

//	typedef struct {
//		int a, b, q, r;         // match score, mismatch penalty and gap open/extension penalty. A gap of size k costs q+k*r
//		int pen_unpaired;       // phred-scaled penalty for unpaired reads
//		int pen_clip;           // clipping penalty. This score is not deducted from the DP score.
//		int w;                  // band width
//		int zdrop;              // Z-dropoff
//
//		int T;                  // output score threshold; only affecting output
//		int flag;               // see MEM_F_* macros
//		int min_seed_len;       // minimum seed length
//		float split_factor;     // split into a seed if MEM is longer than min_seed_len*split_factor
//		int split_width;        // split into a seed if its occurence is smaller than this value
//		int max_occ;            // skip a seed if its occurence is larger than this value
//		int max_chain_gap;      // do not chain seed if it is max_chain_gap-bp away from the closest seed
//		int n_threads;          // number of threads
//		int chunk_size;         // process chunk_size-bp sequences in a batch
//		float mask_level;       // regard a hit as redundant if the overlap with another better hit is over mask_level times the min length of the two hits
//		float chain_drop_ratio; // drop a chain if its seed coverage is below chain_drop_ratio times the seed coverage of a better chain overlapping with the small chain
//		int max_ins;            // when estimating insert size distribution, skip pairs with insert longer than this value
//		int max_matesw;         // perform maximally max_matesw rounds of mate-SW for each end
//		int8_t mat[25];         // scoring matrix; mat[0] == 0 if unset
//	} mem_opt_t;
	/// <summary>
	/// The options for bwa mem, note not all are currently useful or meaningful with this current class, which
	/// allows the user to resolve paired end read mis-matches.
	/// 
	/// Danger: Changes to the scoring matrix likely will force the user to call bwa_fill_scmat, and so this needs to be wrapped better.
	/// </summary>
	[StructLayout(LayoutKind.Sequential)]
	//TODO: The scoring matix likely causes blittable types, and so 
	internal struct mem_opt_t {
		/// <summary>
		/// Match score, default 1
		/// </summary>
		int a;
		/// <summary>
		/// Mismatch penalty, default 4
		/// </summary>
		int b;
		/// <summary>
		/// Deletion open penalty. default 6, A gap of size k costs q+k*r
		/// </summary>
		int o_del;
		/// <summary>
		/// Gap extension penalty. defualt 1,  A gap of size k costs q+k*r
		/// </summary>
		int e_del;
        /// <summary>
        /// Deletion open penalty. default 6, A gap of size k costs q+k*r
        /// </summary>
        int o_ins;
        /// <summary>
        /// Gap extension penalty. defualt 1,  A gap of size k costs q+k*r
        /// </summary>
        int e_ins;


        /// <summary>
		/// phred-scaled penalty for unpaired reads, default 17
		/// </summary>
		int pen_unpaired;       // phred-scaled penalty for unpaired reads
		/// <summary>
		/// clipping penalty. default 5, This score is not deducted from the DP score.
		/// </summary>
		int pen_clip5;
        int pen_clip3;
		/// <summary>
		/// band width, defualt 100
		/// </summary>
		int w;                  
		/// <summary>
		/// The Z-dropoff. Default 100
		/// </summary>
		int zdrop;    

        ulong max_mem_intv;

		/// <summary>
		/// Output score threshold; only affecting output. Default 30.
		/// </summary>
		int T;                  // output score threshold; only affecting output
		/// <summary>
		/// The flag, see MEM_F_* macros, Default 0.
		/// </summary>
		int flag;               // see MEM_F_* macros
		/// <summary>
		/// The minimum seed length.  Default 19.
		/// </summary>
		int min_seed_len;       // minimum seed length

        int min_chain_weight;
        int max_chain_extend;

		/// <summary>
		/// The split_factor. Default 1.5. Split into a seed if MEM is longer than min_seed_len*split_factor
		/// </summary>
		float split_factor;     // split into a seed if MEM is longer than min_seed_len*split_factor
		/// <summary>
		/// The split_width. Default 10 Split into a seed if its occurence is smaller than this value
		/// </summary>
		int split_width;        // split into a seed if its occurence is smaller than this value
		/// <summary>
		/// Default 10000
		/// </summary>
		int max_occ;            // skip a seed if its occurence is larger than this value
		/// <summary>
		/// The max_chain_gap. Default 10000
		/// </summary>
		int max_chain_gap;      // do not chain seed if it is max_chain_gap-bp away from the closest seed
		/// <summary>
		/// The number of threads. Default 1.
		/// </summary>
		int n_threads;          // number of threads
		/// <summary>
		/// Process chunk_size-bp sequences in a batch.  Default 10000000
		/// </summary>
		int chunk_size;         // process chunk_size-bp sequences in a batch
		/// <summary>
		/// The mask_level. Default 0.50, regard a hit as redundant if the overlap with another better hit is over mask_level times the min length of the two hits
		/// </summary>
		float mask_level;       // regard a hit as redundant if the overlap with another better hit is over mask_level times the min length of the two hits
		/// <summary>
		/// The chain_drop_ratio. Default 0.50. Drop a chain if its seed coverage is below chain_drop_ratio times the seed coverage of a better chain overlapping with the small chain
		/// </summary>
		float drop_ratio; // drop a chain if its seed coverage is below chain_drop_ratio times the seed coverage of a better chain overlapping with the small chain
	
        float XA_drop_ratio;
        float mask_level_redun;
        float mapQ_coef_len;
        /// <summary>
		/// Default - 10000 when estimating insert size distribution, skip pairs with insert longer than this value
		/// </summary>
		int max_ins;            // when estimating insert size distribution, skip pairs with insert longer than this value
		/// <summary>
		/// Perform maximally max_matesw rounds of mate-SW for each end, defualt 100
		/// </summary>
		int max_matesw;         // perform maximally max_matesw rounds of mate-SW for each end	
		//int8_t mat[25];         // scoring matrix; mat[0] == 0 if unset

        int max_XA_hits, max_XA_hits_alt;
		/// <summary>
		/// The mat.
		/// </summary>
		[MarshalAs(UnmanagedType.ByValArray, SizeConst=25)] 
		public byte[] mat;

	}





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

        internal int is_shm;
        internal long l_mem;
        internal IntPtr mem;
	} 

//	typedef struct {
//		int64_t offset;
//		int32_t len;
//		int32_t n_ambs;
//		uint32_t gi;
//		char *name, *anno;
//	} bntann1_t;

	[StructLayout(LayoutKind.Sequential)]
	internal struct bntann1_t{
		internal long offset;
		internal int len;
		internal int n_ambs;
		internal uint gi;
        internal int ia_alt;
		/// <summary>
		/// char*
		/// </summary>
		internal IntPtr name;
		/// <summary>
		/// char*
		/// </summary>
		internal IntPtr anno;

	}


	//			//typedef struct {
//				int64_t l_pac;
//				int32_t n_seqs;
//				uint32_t seed;
//				bntann1_t *anns; // n_seqs elements
//				int32_t n_holes;
//				bntamb1_t *ambs; // n_holes elements
//				FILE *fp_pac;
	//		} bntseq_t;
	[StructLayout(LayoutKind.Sequential)]
	internal struct bntseq_t{
		internal long l_pac;
		internal int n_seqs;
		internal uint seed;
		internal IntPtr anns; // n_seqs elements
		internal int n_holes;
		internal IntPtr ambs; // n_holes elements
		internal IntPtr fp_pac;
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

        internal int rid;
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

        internal int alt_sc;

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
		/// Index of the parent hit shadowing the current hit; <0 if primary
		/// </summary>
		internal int secondary;  // index of the parent hit shadowing the current hit; <0 if primary

        internal int secondary_all;
        /// <summary>
        /// Length of the starting seed.
        /// </summary>
        internal int sedlen0;

       
        internal int n_compAndis_alt;
        internal float frac_rep;
        internal ulong hash;
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
		/// is_rev:1, isAlt:1, mapq:8, NM:22;
		/// </summary>
		internal uint isRevAndisAltAndMapQAndNM;
		/// <summary>
		/// number of cigar operations
		/// </summary>
		internal int n_cigar;
		/// <summary>
		/// Pointer to cigar in bam encoding, this is a uint32
		/// </summary>
		internal IntPtr cigar;

        internal IntPtr XA;

		internal int score;
		internal int sub;
        internal int alt_sc;
	}


	//	#define __KSEQ_TYPE(type_t)						\
//typedef struct {							\
//kstring_t name, comment, seq, qual;		\
//int last_char;							\
//kstream_t *f;							\
//} kseq_t;



}

