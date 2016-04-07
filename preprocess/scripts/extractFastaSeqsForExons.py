#!/usr/bin/env python

import sys
import os
import string

# reads from exported environment variable
ExTraMapperPath=os.environ['EXTRAMAPPER_DIR']
sys.path.append(ExTraMapperPath+"/scripts")
from ensemblUtils import parse_organism_GTF
from ensemblUtils import extract_fasta_files_for_exons

#   Testing functionalities
def main(argv):
	indir=argv[1]
	outdir=argv[2]
	orgId1="org1"
	orgId2="org2"

	if not os.path.exists(outdir): os.makedirs(outdir)

	infilename=indir+"/ensemblDownloads/org1.gtf.gz"
#	outfilename=orgId1+"_exons.bed"
#	exonDic1=extract_exon_coordinates_bedFile(orgId1,infilename,outfilename)
	geneDic1,transcriptDic1,exonDic1,infoDic1=parse_organism_GTF(orgId1, infilename, "None")

	infilename=indir+"/ensemblDownloads/org2.gtf.gz"
#	outfilename=orgId2+"_exons.bed"
#	exonDic2=extract_exon_coordinates_bedFile(orgId2,infilename,outfilename)
	geneDic2,transcriptDic2,exonDic2,infoDic2=parse_organism_GTF(orgId2, infilename, "None")
	#

	#
	refGD=indir+"/genomedataArchives/org1" # reference genomedata archive
	typ="allExon"
	fivePrimeFlank=0; threePrimeFlank=0
	outfilename=outdir+"/"+orgId1+"_"+typ+".fasta"
	extract_fasta_files_for_exons(refGD,exonDic1,typ,fivePrimeFlank,threePrimeFlank,outfilename)
	typ="allExonPlusMinus"
	fivePrimeFlank=10; threePrimeFlank=10
	outfilename=outdir+"/"+orgId1+"_"+typ+"_"+str(fivePrimeFlank)+"-5p"+"_"+str(threePrimeFlank)+"-3p"+".fasta"
	extract_fasta_files_for_exons(refGD,exonDic1,typ,fivePrimeFlank,threePrimeFlank,outfilename)
	typ="exonIntron"
	outfilename=outdir+"/"+orgId1+"_"+typ+"_"+str(fivePrimeFlank)+"-5p"+"_"+str(threePrimeFlank)+"-3p"+".fasta"
	extract_fasta_files_for_exons(refGD,exonDic1,typ,fivePrimeFlank,threePrimeFlank,outfilename)
	typ="intronExon"
	outfilename=outdir+"/"+orgId1+"_"+typ+"_"+str(fivePrimeFlank)+"-5p"+"_"+str(threePrimeFlank)+"-3p"+".fasta"
	extract_fasta_files_for_exons(refGD,exonDic1,typ,fivePrimeFlank,threePrimeFlank,outfilename)
	
	##

	refGD=indir+"/genomedataArchives/org2" # reference genomedata archive
	typ="allExon"
	fivePrimeFlank=0; threePrimeFlank=0
	outfilename=outdir+"/"+orgId2+"_"+typ+".fasta"
	extract_fasta_files_for_exons(refGD,exonDic2,typ,fivePrimeFlank,threePrimeFlank,outfilename)
	typ="allExonPlusMinus"
	fivePrimeFlank=10; threePrimeFlank=10
	outfilename=outdir+"/"+orgId2+"_"+typ+"_"+str(fivePrimeFlank)+"-5p"+"_"+str(threePrimeFlank)+"-3p"+".fasta"
	extract_fasta_files_for_exons(refGD,exonDic2,typ,fivePrimeFlank,threePrimeFlank,outfilename)
	typ="exonIntron"
	outfilename=outdir+"/"+orgId2+"_"+typ+"_"+str(fivePrimeFlank)+"-5p"+"_"+str(threePrimeFlank)+"-3p"+".fasta"
	extract_fasta_files_for_exons(refGD,exonDic2,typ,fivePrimeFlank,threePrimeFlank,outfilename)
	typ="intronExon"
	outfilename=outdir+"/"+orgId2+"_"+typ+"_"+str(fivePrimeFlank)+"-5p"+"_"+str(threePrimeFlank)+"-3p"+".fasta"
	extract_fasta_files_for_exons(refGD,exonDic2,typ,fivePrimeFlank,threePrimeFlank,outfilename)

	return

if __name__ == "__main__":
	main(sys.argv)
#

