###Introduction

This library is a C# wrapper around the bwa mem program API.  This code allows you to call BWA from C# (or python or F#) and obtain fully typed SAMAlignment object compatible with .NET Bio classes.  It differs from previous wrappers as it does not require the user to run the index command prior to calling this program, and also does not to the interop through a C++ layer.  Bwa mem was originally written by Heng Li and is available at: http://bio-bwa.sourceforge.net/.  BWA was released under GPLv3, which I trust was a good choice, and so this license also applies to this code. 

This C# library comes as a .dll file.  This file depends on a shared library that is compiled by gcc as 64 bit.  To use on 32 bit systems the user must alter the make file appropriately.

Because the original BWA program has a dependencies on linux, this library only works on posix systems that are little endian and are executing as 64 bit programs.  It therefore must also rely on the gnu build tool chain, so cannot be expected to run as smoothly as CLS compliant code.  Installation issues may occur.

###Getting started

Example coming soon...

###Citing the original BWA

* Li H. and Durbin R. (2009) Fast and accurate short read alignment with
 Burrows-Wheeler transform. *Bioinformatics*, **25**, 1754-1760. [PMID:
 [19451168][1]]. (if you use the BWA-backtrack algorithm)

* Li H. and Durbin R. (2010) Fast and accurate long-read alignment with
 Burrows-Wheeler transform. *Bioinformatics*, **26**, 589-595. [PMID:
 [20080505][2]]. (if you use the BWA-SW algorithm)

* Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs
 with BWA-MEM. [arXiv:1303.3997v2][3] [q-bio.GN]. (if you use the BWA-MEM
 algorithm or the **fastmap** command, or want to cite the whole BWA package)

Please note that the last reference is a preprint hosted at [arXiv.org][4]. I
do not have plan to submit it to a peer-reviewed journal in the near future.

[1]: http://www.ncbi.nlm.nih.gov/pubmed/19451168
[2]: http://www.ncbi.nlm.nih.gov/pubmed/20080505
[3]: http://arxiv.org/abs/1303.3997
[4]: http://arxiv.org/
