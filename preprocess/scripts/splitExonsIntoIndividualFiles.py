#!/usr/bin/env python
import sys

def main(argv):
	infilename=argv[1]
	outdir=argv[2]
	whichCol=int(argv[3])-1
	fileSuffix=argv[4]
	infile=open(infilename,'r')
	lastExon="dummy"
	outfile=open("dummy.txt",'w')
        for line in infile:
                newExon=line.rstrip().split()[whichCol] # where exon name is
		if newExon!=lastExon:
			outfile.close()
			outfile=open(outdir+"/"+newExon+fileSuffix,'w')
		#
		outfile.write(line)
		lastExon=newExon
	#
	outfile.close()

	return

if __name__ == "__main__":
        main(sys.argv)
#

