#!/usr/bin/env python
import subprocess,sys,struct,os,time,shutil



if len(sys.argv) < 4:
	sys.exit( 'Usage : [Read file] [k] [nks (min occurence of kmer)]' )

debut=time.time()

readFile=sys.argv[1]
k=sys.argv[2]
nks=sys.argv[3]
output='out'


print "PHASE 1 : Get Kmers"
subprocess.call(['dsk','-file', readFile ,'-kmer-size', k,'-abundance-min',nks,'-max-memory','5000','-out','reads'])
subprocess.call(['dsk2ascii','-file', 'reads.h5','-out', 'kmers'])


print "PHASE 2 : Get Unitigs"
subprocess.call(['./bcalm', 'kmers','unitig.dot','6'])
#shutil.rmtree('.bcalmtmp')


print "PHASE 3 : Map Reads on graph"
subprocess.call(['/usr/bin/time','-f',"\'%e %P %M\'",'./bgreat','-r',readFile,'-g','unitig.dot','-k',k,'-t','4'])
#os.remove('paths')

print "PHASE 4 : Map Reads on big unitigs with bowtie2"
#see sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/
subprocess.call(['./getLargeUnitigs','unitig.dot','big.fa','100'])
subprocess.call(['./bowtie2-build','big.fa','index','-q'])
subprocess.call(['./bowtie2','-r','--very-fast','-x','index','-U','notAligned.fa','-t','-p','4','-S','out.sam'])

os.remove('paths')
#os.remove('big.fa')
os.remove('kmers')
os.remove('largest_bucket.dot')
#os.remove('noOverlap.fa')
#os.remove('notAligned.fa')
os.remove('out.sam')
os.remove('reads.h5')
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
