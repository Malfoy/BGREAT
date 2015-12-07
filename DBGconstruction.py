#!/usr/bin/env python
import subprocess,sys,struct,os,time,shutil



if len(sys.argv) < 5:
	sys.exit( 'Usage : [Read file] [k] [nks (min occurence of kmer)] [output file]' )

debut=time.time()

readFile=sys.argv[1]
k=sys.argv[2]
nks=sys.argv[3]
output=sys.argv[4]


print "PHASE 1 : Get Kmers"
subprocess.call(['dsk','-file', readFile ,'-kmer-size', k,'-abundance-min',nks,'-max-memory','5000','-out','reads'])
subprocess.call(['dsk2ascii','-file', 'reads.h5','-out', 'kmers'])


print "PHASE 2 : Get Unitigs"
subprocess.call(['./bcalm', 'kmers',output,'6'])

print 'Cleaning'
shutil.rmtree('.bcalmtmp')
os.remove('kmers')
os.remove('reads.h5')

fin=time.time();
print fin-debut
print "secondes"
