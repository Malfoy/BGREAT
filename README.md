Bgreat
=====
De Bruijn Graph Read Mapping Tool

Compilation
====
	$make

Usage
====
	$bgreat  reads_file k unitigs_file number_missmatchs_allowed number_threads path_file no_overlap_file not_aligned_file

The default value are:

k:31

unitigs_file: unitig.dot

number_missmatchs_allowed: 2

number_threads: 1

path_file: paths

no_overlap_file: noOverlap.fa

not_aligned_file: notAligned.fa


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

0.1.2.

1.1.3.4.-6.

The first integer is the starting position of the read in the first unitig.

The other integers are the unitig numbers in the order they appear in the read, the minus sign mean that the reverse complement of the unitig should be used.



!!!Warning : this version is linux intended and still in testing!!!
