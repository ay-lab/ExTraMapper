#!/usr/bin/env python
##############################################################################
### To use the functions in this lib simply import this python module using
### import ensemblUtils
### Then you'll able able to call functions with the proper arguments using
### returnVal=ensemblUtils.func1(arg1,arg2)
##############################################################################
##############################################################################
import sys
import os
import string
import math
import gzip
import cPickle as pickle

# reads from exported environment variable
ExTraMapperPath=os.environ['EXTRAMAPPER_DIR']
sys.path.append(ExTraMapperPath+"/scripts")
from ensemblUtils import *

#   Testing functionalities
def main(argv):
	indir=argv[1]
	orgId1="org1"; orgId2="org2";
	refGD1=indir+"/genomedataArchives/org1"
	refGD2=indir+"/genomedataArchives/org2"


#	outdir="GTFsummaries"; 
	if len(argv)==2:
		return

	outdir=argv[2]
	os.system("mkdir -p "+outdir)
		
	infilename=indir+"/ensemblDownloads/org1.gtf.gz"
	geneDic1,transcriptDic1,exonDic1,infoDic1=parse_organism_GTF(orgId1, infilename, outdir)
	
	infilename=indir+"/ensemblDownloads/org2.gtf.gz"
	geneDic2,transcriptDic2,exonDic2,infoDic2=parse_organism_GTF(orgId2, infilename, outdir)

	## these two files were downloaded by hand selecting columns from Ensembl's Biomart
	## I weren't able to redo the same column selections recently so I decided to switch to
	## parsing the orthology information from readily available Ensembl files like below ones:
	## ftp://ftp.ensembl.org/pub/release-80/mysql/ensembl_mart_80/
	## hsapiens_gene_ensembl__homolog_mmus__dm.txt.gz
	#infilename="/projects/b1017/shared/Ensembl-files/Ensembl-human-GRCh38-to-mouse-GRCm38.p3.txt.gz"
	#genePairsHumanToMouse=parse_ensembl_gene_pairings(infilename)
	#infilename="/projects/b1017/shared/Ensembl-files/Ensembl-mouse-GRCm38.p3-to-human-GRCh38.txt.gz"
	#genePairsMouseToHuman=parse_ensembl_gene_pairings(infilename)
	#consistency_check(genePairsHumanToMouse,genePairsMouseToHuman)
	## if consistency check is ok then just use one side. This is OK for one2one mappings.
	#genePairsDic=genePairsHumanToMouse
	#pickle_one2one_genePairs_allInfo(genePairsDic,geneDic1,geneDic2,exonDic1,exonDic2,transcriptDic1,transcriptDic2,outdir) 
	
	infilename=indir+"/ensemblDownloads/org1_homolog_org2.txt.gz"
	proteinToGeneDic,genePairsDic,proteinPairsDic=parse_ensembl_geneAndProtein_pairings(infilename,{},{})
	print ["1",len(proteinToGeneDic),len(genePairsDic),len(proteinPairsDic)]

	infilename=indir+"/ensemblDownloads/org2_homolog_org1.txt.gz"
	proteinToGeneDic,genePairsDic,proteinPairsDic=parse_ensembl_geneAndProtein_pairings(infilename,proteinToGeneDic,proteinPairsDic)
	print ["2",len(proteinToGeneDic),len(genePairsDic),len(proteinPairsDic)]

	
	exonDic1=assign_firstMidLast_exon_counts(exonDic1,transcriptDic1)
	exonDic2=assign_firstMidLast_exon_counts(exonDic2,transcriptDic2)

	typ="allExonPlusMinus"
	outfilename="None"
	fivePrimeFlank=12; threePrimeFlank=12
	exonDic1=extract_conservation_stats_for_exons(refGD1,exonDic1,typ,fivePrimeFlank,threePrimeFlank,outfilename)
	exonDic2=extract_conservation_stats_for_exons(refGD2,exonDic2,typ,fivePrimeFlank,threePrimeFlank,outfilename)
 
	outdir=argv[2] # overwrite previous summaries
	os.system("mkdir -p "+outdir)
	print_some_summary(orgId1, geneDic1,transcriptDic1,exonDic1,{}, outdir)
	print_some_summary(orgId2, geneDic2,transcriptDic2,exonDic2,{}, outdir)

#	outdir="perGenePairExonLists"
	if len(argv)==3:
		return

	outdir=argv[3]
	os.system("mkdir -p "+outdir)

	outfilename=outdir+"/genePairsSummary-one2one.txt"
	print_one2one_genePairs(genePairsDic,geneDic1,geneDic2,outfilename) # either way is ok since one2one

	geneOnlyOrthoDic1,transcriptOnlyOrthoDic1,exonOnlyOrthoDic1, geneOnlyOrthoDic2,transcriptOnlyOrthoDic2,exonOnlyOrthoDic2=pickle_one2one_genePairs_allInfo(genePairsDic,geneDic1,geneDic2,exonDic1,exonDic2,transcriptDic1,transcriptDic2,outdir) 

	print [len(geneDic1), len(geneDic2)]
	print len(geneOnlyOrthoDic1)
	print len(geneOnlyOrthoDic2)

	outdir=argv[2]+"/onlyOrthologAndCodingGenes"
	os.system("mkdir -p "+outdir)
	print outdir
	print_some_summary(orgId1, geneOnlyOrthoDic1,transcriptOnlyOrthoDic1,exonOnlyOrthoDic1,{}, outdir)
	print_some_summary(orgId2, geneOnlyOrthoDic2,transcriptOnlyOrthoDic2,exonOnlyOrthoDic2,{}, outdir)

	#
#	print_one2one_exonListPairs(genePairsDic,geneDic1,geneDic2,exonDic1,exonDic2,orgId1,orgId2,outdir)
#	print_one2one_transcriptListPairs(genePairsDic,geneDic1,geneDic2,transcriptDic1,transcriptDic2,orgId1,orgId2,outdir)

	return

if __name__ == "__main__":
	main(sys.argv)

