Bgreat
=====
De Bruijn Graph Read Mapping Tool

Compilation
====
	$make

Usage
====
	$bgreat -r read_files -k k_value -g unitig_file -m n_missmatch -t n_thread -p path_file -o no_overlap_file
	-a not_aligned_file

You can give multiple fasta file in input by separating the file names with a coma

	$./bgreat -r reads.fa,reads2.fa -k 4 -g unitig.fa

If you use fastq files (with no splited line) use the -q option

The default value are:


The k value of your De Bruijn graph
"k_value":31


The file containing your unitigs
"unitig_file": unitig.fa


The number of missmatches allowed in the alignment
"n_missmatch_allowed": 2


The number of thread used for mapping
"n_thread": 1


The paths ouput file
"path_file": paths


The output file with read without overlap, that may be mapped by a regular aligner on large unitigs
"no_overlap_file": noOverlap.fa


The outputfile with read that have an overlap but could not be mapped by BGREAT
"not_aligned_file": notAligned.fa


Graph computation
====


We propose the following method to compute the unitigs from a set of reads:

Unitig computation with BCALM2 (see https://github.com/GATB/bcalm) :

	&./bcalm -in reads.fa -k 31 -abundance 3
	&./bglue -in unitigs.h5 -k 31

Output
====
You can try

	$./bgreat -r reads.fa -k 4 -g unitig.fa

The Output file will look like something like this

>1miss

0.1.2.

>nomiss

0.3.4.-6.


The first integer is the starting position of the read in the first unitig.

The other integers are the unitig numbers in the order they appear in the read, the minus sign mean that the reverse complement of the unitig should be used.
