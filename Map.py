#!/usr/bin/env python
import subprocess,sys,struct,os,time,shutil



if len(sys.argv) < 4:
	sys.exit( 'Usage : [Read file] [k] [nks (min occurence of kmer)]' )

debut=time.time()

readFile=sys.argv[1]
k=sys.argv[2]
nks=sys.argv[3]
output='out'
readsize='100'
errorRate='0.01'
unitigFile=sys.argv[4]

#print "PHASE 0 : Get Reads"
#subprocess.call(['./mutareads', 'ecoliref.fasta', 'reads', nbReads , readsize, errorRate, '0', '0','-snp','0.001'])
#subprocess.call(['./ReadSim','-file','CELEGANSREF.fa','-out','reads.fa','-nb-read', nbReads ,'-min-read-size', readsize,'-error', errorRate])


#print "PHASE 1 : Get Kmers"
#subprocess.call(['dsk','-file', readFile ,'-kmer-size', k,'-abundance-min',nks,'-max-memory','10000','-out','reads'])
#subprocess.call(['dsk2ascii','-file', 'reads.h5','-out', 'kmers'])
#
#
#print "PHASE 2 : Get Unitigs"
#subprocess.call(['bcalm', 'kmers','unitig.dot','6'])
#shutil.rmtree('.bcalmtmp')


print "PHASE 3 : Map Reads on graph"
subprocess.call(['./BGREAT',k,readFile,unitigFile,'2','1'])
#os.remove('paths')

print "PHASE 4 : Map Reads on big unitigs"
subprocess.call(['./getBigUnitig',unitigFile,'big.fa',readsize])
subprocess.call(['bowtie2-build','big.fa','index','-q'])
subprocess.call(['bowtie2','-f','--very-fast','-x','index','-U','notAligned.fa','-t','-p','4','-S','out.sam'])

os.remove('paths')
os.remove('big.fa')
#os.remove('kmers')
#os.remove('largest_bucket.dot')
os.remove('noOverlap.fa')
os.remove('notAligned.fa')
os.remove('out.sam')
#os.remove('reads.h5')
#os.remove('unitig.dot')
os.remove('index.1.bt2')
os.remove('index.2.bt2')
os.remove('index.3.bt2')
os.remove('index.4.bt2')
os.remove('index.rev.1.bt2')
os.remove('index.rev.2.bt2')




print "The End"
fin=time.time();
print fin-debut
print "secondes"
