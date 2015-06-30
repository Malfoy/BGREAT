Bgreat
=====
De Bruijn Graph Read Mapping Tool

Compilation
====
	$make

Usage
====
	$bgreat  read_file k unitig_file n_missmatch n_thread path_file no_overlap_file not_aligned_file

 The default value are:
 

The k value of your De Bruijn graph

"k":31
 

The file containing your unitigs

"unitig_file": unitig.dot
 

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


Workflow
====


We propose the following workflow to compute unitigs from a set of reads:


Kmer counting with DSK (see http://minia.genouest.org/dsk/) :

	$dsk -file reads.fasta -kmer-size k -abundance-min c -out reads

	$dsk2ascii -file reads.h5 -out kmers


Unitig computation with BCALM (see https://github.com/Malfoy/bcalm) :

	$bcalm kmers unitigs.dot 6

Output
====
You can try

	$./bgreat reads.fa 4 unitig.dot

The Output file will look like something like this

0.1.2.5.

1.1.3.4.-6.0.

The first integer is the starting position of the read in the first unitig.

The last integer is the number of char of ther last unitig that is used for the alignment, 0 mean that the whole unitig is used.

The other integers are the unitig numbers in the order they appear in the read, the minus sign mean that the reverse complement of the unitig should be used.




!!!Warning : this version is linux intended and still in testing!!!
