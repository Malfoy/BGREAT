Bgreat : De Bruijn Graph Read Mapping Tool

Compilation :
make

Usage :
bgreat k reads.fasta unitigs.dot number_missmatchs_allowed number_threads

We propose the following worflow to compute unitigs from a set of reads:

Kmer counting with DSK (see http://minia.genouest.org/dsk/):
dsk -file reads.fasta -kmer-size k -abundance-min c -out reads
dsk2ascii -file reads.h5 -out kmers

Unitig computation with BCALM (see github.com/Malfoy/bcalm) :
bcalm kmers unitigs.dot 6 
