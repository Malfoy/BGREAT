#!/usr/bin/env python
import subprocess,sys,struct,os,time,shutil



if len(sys.argv) < 5:
	sys.exit( 'Usage : [Read file] [unitig file] [output file] [number thread]' )

debut=time.time()

readFile=sys.argv[1]
unitigFile=sys.argv[2]
outputFile=sys.argv[3]
nthread=sys.argv[4]



print "PHASE 3 : Map Reads on graph"
subprocess.call(['/usr/bin/time','-f',"\'%e %P %M\'",'./bgreat','-r',readFile,'-g',unitigFile,'-k',k,'-t',nthread])


print "PHASE 4 : Map Reads on big unitigs with bowtie2"
#see sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/
subprocess.call(['./getLargeUnitigs',unitigFile,'big.fa','100'])
subprocess.call(['./bowtie2-build','big.fa','index','-q'])
subprocess.call(['./bowtie2','-r','--very-fast','-x','index','-U','notAligned.fa','-t','-p',nthread,'-S','out.sam'])

os.remove('big.fa')
fin=time.time();
print fin-debut
print "secondes"
