Bgreat
=====
De Bruijn Graph Read Mapping Tool

Compilation
====
  $make

Usage
====
	$bgreat k reads.fasta unitigs.dot number_missmatchs_allowed number_threads

We propose the following workflow to compute unitigs from a set of reads:


Kmer counting with DSK (see http://minia.genouest.org/dsk/):

	$dsk -file reads.fasta -kmer-size k -abundance-min c -out reads

	$dsk2ascii -file reads.h5 -out kmers


Unitig computation with BCALM (see github.com/Malfoy/bcalm) :

	$bcalm kmers unitigs.dot 6

Output
====
You can try 

	$./bgreat 4 reads.fa unitig.dot  1 1

The Output file will look like something like this

0.1.2.

1.1.3.4.-6.

The first integer is the starting position of the read in the first Unitig.
The other integer are the unitig number in the order they appears in the read, the minus sign mean that the reverseComplement of the unitig should be used.



!!!Warning : this version is linux intended and still in testing!!!
