#!/usr/bin/env python
import sys
import gzip

def main(argv):
	infilename=argv[1]
	outdir=argv[2]
	lastCh="dummy"
	outfile=open("dummy.txt",'w')
	if infilename.endswith(".gz"):
		infile=gzip.open(infilename,'r')
	else:
		infile=open(infilename,'r')
	
	for l in infile:
		if l.startswith(">"):
			chrStr=l.split()[0]
			ch=chrStr[1:]
			outfile.close()
			outfile=open(outdir+"/"+ch+".fa",'w')
			outfile.write(">"+ch+"\n")
		else:
			outfile.write(l)
	#
	infile.close()
	outfile.close()
	
	return

if __name__ == "__main__":
        main(sys.argv)
#

