using System;


namespace Bio.BWA.MEM
{
	class MainClass
	{
		static void Main (string[] args)
		{
			Console.WriteLine ("Hello World");
			BWA bwa = new BWA ("/Users/ndelaney/bwa-0.7.5a/TestData/MT.fasta");
			Sequence s = new Sequence (DnaAlphabet.Instance, "CGCATTCCTACTACTCAACTTAAACTCCAGCACCACGACCCTACTACTATCTCGCACCTGCTTTT");
			bwa.AlignSequence (s);


		}
	}
}	

