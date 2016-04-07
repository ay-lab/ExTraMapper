#!/usr/bin/env python

import sys
import os
import string

# reads from exported environment variable
ExTraMapperPath=os.environ['EXTRAMAPPER_DIR']
sys.path.append(ExTraMapperPath+"/scripts")
from ensemblUtils import parse_organism_GTF
#from ensemblUtils import EnsemblExon


def main(argv):
	"""
	This script will generate some output of the form below which is required by Exalign: 

	#refGene.name	refGene.chrom	refGene.strand	refGene.cdsStart	refGene.cdsEnd	refGene.exonCount	refGene.exonStarts	refGene.exonEnds	refGene.name2	refLink.protAc

	ENST00000618887	chr17	-	28348120	28357448	11	28347661,28347661,28349083,28350438,28350766,28351664,28352912,28353241,28353695,28354488,28355795,28357288,	28348231,28348231,28349162,28350563,28350792,28351800,28353019,28353316,28353791,28354531,28355876,28357463,	ENSG00000004142	ENSP00000477665

	ENST00000618164	chr11	+	86306412	86345635	5	86306333,86306416,86337379,86344603,86345584,	86306414,86306482,86337530,86344721,86345662,	ENSG00000149196	ENSP00000482151

	"""

	indir=argv[1]
	outdir=argv[2]
	orgId1="org1"
	orgId2="org2"

	if not os.path.exists(outdir): os.makedirs(outdir)

	infilename=indir+"/ensemblDownloads/org1.gtf.gz"
	geneDic1,transcriptDic1,exonDic1,infoDic1=parse_organism_GTF(orgId1, infilename, "None")

	infilename=indir+"/ensemblDownloads/org2.gtf.gz"
	geneDic2,transcriptDic2,exonDic2,infoDic2=parse_organism_GTF(orgId2, infilename, "None")

	outfile1=open(outdir+"/org1.rf",'w')
	outfile2=open(outdir+"/org2.rf",'w')
	outfiles=[outfile1, outfile2]
	dics=[transcriptDic1,transcriptDic2]
	for x in [0,1]:
		transcriptDic=dics[x]
		outfile=outfiles[x]
		for t in transcriptDic:
			tEntry=transcriptDic[t]
			if tEntry.basicInfoDic["transcript_biotype"]!="protein_coding":
				continue
			id=t; 
			ch="chr"+tEntry.basicInfoDic["chromosome"]
			strand=tEntry.basicInfoDic["strand"]
			codingS=-1
			codingE=-1
			noExons=len(tEntry.exons)
			exonStarts=",".join([str(e[0]) for e in tEntry.exons])+","
			exonEnds=",".join([str(e[1]) for e in tEntry.exons])+","
			geneName=tEntry.basicInfoDic["gene_id"]
			proteinName=tEntry.protein_id

			if len(tEntry.start_codon)==2:
				codingS=tEntry.start_codon[0]
			if len(tEntry.stop_codon)==2:
				codingE=tEntry.stop_codon[1]
			if codingS==-1 or codingE==-1:
				maxC=0; minC=9999999999
				for tempS,tempE in tEntry.codingExons:
					if tempS==-1 and tempE==-1:
						continue
					elif tempS==-1:
						tempS=tempE
					elif tempE==-1:
						tempE=tempS
					if max(tempS,tempE)>maxC:
						maxC=max(tempS,tempE)
					if min(tempS,tempE)<minC:
						minC=min(tempS,tempE)
				#
				codingS=minC
				codingE=maxC
			#
			if strand=="-":
				exonStarts=",".join([str(e[0]) for e in list(reversed(tEntry.exons))])+","
				exonEnds=",".join([str(e[1]) for e in list(reversed(tEntry.exons))])+","

			outfile.write("%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n" % (id,ch,strand,codingS,codingE,noExons,exonStarts,exonEnds,geneName,proteinName))
		#
		outfile.close()
	#
#
#		summary="chr"+self.basicInfoDic["chromosome"]+"\t"+str(self.basicInfoDic["start_coord"])+"\t"+\
#			str(self.basicInfoDic["end_coord"])+"\t"+self.basicInfoDic["strand"]+"\t"+self.basicInfoDic["transcript_id"]+"\t"+\
#			self.basicInfoDic["transcript_name"]+"\t"+self.basicInfoDic["transcript_biotype"]+"\t"+\
#			self.basicInfoDic["gene_id"]+"\t"+self.basicInfoDic["gene_name"]+"\t"+\
#			",".join([e[2] for e in self.exons])+"\t"+",".join(self.exon_types)+"\t"+\
#			",".join([str(e[0]) for e in self.exons])+"\t"+",".join([str(e[1]) for e in self.exons])+"\t"+\
#			",".join([str(e[0]) for e in self.codingExons])+"\t"+",".join([str(e[1]) for e in self.codingExons])+"\t"+\
#			",".join([str(i) for i in self.start_codon])+"\t"+",".join([str(i) for i in self.stop_codon])+"\t"+\
#			",".join([str(e[0]) for e in self.UTRs])+"\t"+",".join([str(e[1]) for e in self.UTRs])+"\t"+\
#			self.protein_id+"\t"+"transcript"

	return
	
if __name__ == "__main__":
	main(sys.argv)
#

