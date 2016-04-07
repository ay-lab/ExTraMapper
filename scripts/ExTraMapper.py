#!/usr/bin/env python

import sys
import os
import string
import operator
import cPickle as pickle
complement = string.maketrans('atcgn', 'tagcn')

# reads from exported environment variable
ExTraMapperPath=os.environ['EXTRAMAPPER_DIR']
sys.path.append(ExTraMapperPath+"/scripts")
from ensemblUtils import *

#import matplotlib
#matplotlib.use('Agg')
#from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
#import matplotlib.pyplot as plt
import numpy as np
#

#### matplotlib fontsize settings
#plt.rcParams['font.size']=17
#plt.rcParams['axes.labelsize']='large'
#plt.rcParams['xtick.labelsize']='medium'
#plt.rcParams['ytick.labelsize']='medium'
#plt.rcParams['figure.subplot.hspace']=0.4
#plt.rcParams['figure.subplot.bottom']=0.12
#plt.rcParams['figure.subplot.left']=0.15
#plt.rcParams['figure.subplot.right']=0.94
#plt.rcParams['figure.subplot.top']=0.92
#########################


def compute_overlap_score_for_exons(exonLength1, exonLength2, liftedOverLength1, overlap):
	r=1 # if exonLength1==liftedOverLength1, meaning liftOver did not change the length
#	print [overlap, liftedOverLength1, exonLength2, exonLength1]
	overlapScore=(2.0*overlap)/(liftedOverLength1+exonLength2)
	if exonLength1<liftedOverLength1:
		r=(1.0*exonLength1)/liftedOverLength1
	elif exonLength1>liftedOverLength1:
		r=(1.0*liftedOverLength1)/exonLength1
	return min(1.0, r*overlapScore)

def parse_liftOver_nonintersecting_file_perExon(infilename,oneToTwo):
	"""
	Parse the .nonintersecting files and determine splice junction losses.
	An example line:  
		chr18	9495767	9495820	+	ENSE00000327880	2	flank0-minMatch0.9-multiples
		chr6	119893048	119893101	+	ENSE00000327880	3	flank0-minMatch0.9-multiples
	"""
	infile=open(infilename,'r')

	if oneToTwo==True:
		org="org2"
		exonDicToUse=exonDic1
		geneEntry=geneDic2[g2]
	else:
		org="org1"
		exonDicToUse=exonDic2
		geneEntry=geneDic1[g1]
	#

	with Genome(genomedataDir+"/"+org) as genome:
		# move to above
		for line in infile:
			words=line.rstrip().split()
			# if exon names are coming from org1, lifted over coords a from org2
			ch,s,e,strand,ex,multiplicity,suffixStr=words
			minMatchStr=suffixStr.split("-")[1][8:] # parse the suffixStr assuming a specific format 
			# don't use the entries coming from partCoding regions as they will not have the splice sites
			if suffixStr.endswith("partCoding"):
				continue
			
			gCh,gS,gE="chr"+str(geneEntry.basicInfoDic["chromosome"]),geneEntry.basicInfoDic["start_coord"],geneEntry.basicInfoDic["end_coord"]
			# discard liftover entries that have nothing to do with the corresponding ortholog gene
			if gCh!=ch:
			#	print ["## nonintersectingNull ###########",ex,ch,s,e,strand,gCh,gS,gE,geneEntry.basicInfoDic["gene_name"]]
				continue
			elif max(0, min(int(e),int(gE))-max(int(s),int(gS)))<1:
			#	print ["## nonintersectingNull2 ###########",ex,ch,s,e,strand,gCh,gS,gE,geneEntry.basicInfoDic["gene_name"]]
				continue
			
			# information about the original exon before lifted over
			origExonEntry=exonDicToUse[ex]
			acceptor1=origExonEntry.acceptor2bp
			donor1=origExonEntry.donor2bp
			avgConsScore=origExonEntry.avgConsScore
			fml1=origExonEntry.firstMidLast #[f,m,l] counts
			orCh,orS=origExonEntry.basicInfoDic["chromosome"],origExonEntry.basicInfoDic["start_coord"]
			orE,orStrand=origExonEntry.basicInfoDic["end_coord"],origExonEntry.basicInfoDic["strand"]
			origExonLen=abs(orS-orE)
			liftedOverLenForEx=abs(int(s)-int(e))

			## Not needed anymore due to above ch==gCh check
			#if ch not in genome:
			#	continue

			sq=genome[ch].seq[int(s)-3:int(e)+2].tostring().lower().upper()
			if strand=="-":
				sqAcc=sq[-2:].lower().translate(complement)[::-1].upper() # reverse complement
				sqDon=sq[:2].lower().translate(complement)[::-1].upper() # reverse complement
			else:
				sqAcc=sq[:2]
				sqDon=sq[-2:]
			#

			f,m,l=fml1
			junctionProblem=False
			if m>0 or (f>0 and l>0):
				if sqAcc!=acceptor1 or sqDon!=donor1:
					junctionProblem=True
				if junctionProblem==True:
				#if True:
					k=acceptor1+"-"+donor1
					k2=sqAcc+"-"+sqDon
					if k not in bigAcceptorDonorSiteDic:
						bigAcceptorDonorSiteDic[k]={}
					if k2 not in bigAcceptorDonorSiteDic[k]:
						bigAcceptorDonorSiteDic[k][k2]=0
					bigAcceptorDonorSiteDic[k][k2]=bigAcceptorDonorSiteDic[k][k2]+1
					
			elif f>0:
				if sqDon!=donor1:
					junctionProblem=True
				if junctionProblem==True:
				#if True:
					if donor1 not in bigDonorSiteDic:
						bigDonorSiteDic[donor1]={}
					if sqDon not in bigDonorSiteDic[donor1]:
						bigDonorSiteDic[donor1][sqDon]=0
					bigDonorSiteDic[donor1][sqDon]=bigDonorSiteDic[donor1][sqDon]+1
			elif l>0:
				if sqAcc!=acceptor1:
					junctionProblem=True
				if junctionProblem==True:
				#if True:
					if acceptor1 not in bigAcceptorSiteDic:
						bigAcceptorSiteDic[acceptor1]={}
					if sqAcc not in bigAcceptorSiteDic[acceptor1]:
						bigAcceptorSiteDic[acceptor1][sqAcc]=0
					bigAcceptorSiteDic[acceptor1][sqAcc]=bigAcceptorSiteDic[acceptor1][sqAcc]+1
			#
			if ex not in nonintersectingExons:
				nonintersectingExons[ex]={}
			key=ch+"-"+str(s)+"-"+str(e)
			# keep the entry with maximum minMatch value
			if key not in nonintersectingExons[ex]:
				# original acceptor and donor sites, first mid last exon counts, liftedover, 
				nonintersectingExons[ex][key]=[minMatchStr,junctionProblem,ch,s,e,strand,sqAcc,sqDon,liftedOverLenForEx]
			else:
				if float(minMatchStr)>float(nonintersectingExons[ex][key][0]):
					nonintersectingExons[ex][key][0]=minMatchStr
			#
#			print ("PotentialSpliceLoss\t%r\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%d\n" % \
#				(junctionProblem,minMatchStr,origExonEntry.get_summary_string(),ch,int(s),int(e),strand,sqAcc,sqDon,liftedOverLenForEx)),
		#
	#
	infile.close()
	return

# no problem with symmetry, just add all unmapped exons from org1 and org2 into unmappedExons dictionary
def parse_liftOver_unmapped_file_perExon(infilename):
	"""
	Parse the .unmapped files.
	An example line:  
		chrX	100593624	100594035	+	ENSE00001952391	#Split	flank0-minMatch1.0-multiples
	"""
	infile=open(infilename,'r')

	for line in infile:
		words=line.rstrip().split()
		ch,s,e,strand,ex,whyUnmapped,suffixStr=words
		minMatchStr=suffixStr.split("-")[1][8:] # parse the suffixStr assuming a specific format 

		if whyUnmapped not in unmappedExons[ex]:
			unmappedExons[ex][whyUnmapped]=0 # just something small, will be overwritten
		#
		if float(unmappedExons[ex][whyUnmapped])<float(minMatchStr):
			unmappedExons[ex][whyUnmapped]=minMatchStr 
#		print ("UnmappedExon\t%s\t%s\t%d\t%d\t%s\t%s\t%s\n" % (ex,ch,int(s),int(e),strand,whyUnmapped,minMatchStr)),
	#
	infile.close()
	return

def parse_liftOver_mapped_file_perExon(infilename,oneToTwo):
	"""
	Parse the .mapped files.
	An example line: 
		chr5	108461232	108461451	+	ENSMUSE00000692363		
		chr5	108461360	108461451	+	ENSE00001875286	
		91	flank0-minMatch0.95-multiples
	"""

	infile=open(infilename,'r')

	for line in infile:
		words=line.rstrip().split()
		# lifted over coordinates for ex1 and original coords for ex2. For human to mouse, this is 1-mouse, 2-human
		ch2,s2,e2,strand2,ex2,ch1,s1,e1,strand1,ex1,dummy,suffixStr=words #ignore the overlap here, it is 1-bp off
		liftedOverLenForEx=abs(int(s1)-int(e1)) # liftover length is always whatever comes second in the file
		codingOverlapScore=0; allOverlapScore=0; 
		minMatchStr=suffixStr.split("-")[1][8:] # parse the suffixStr assuming a specific format 


		if oneToTwo==True:
			# discard those that has nothing to do with this gene
			if ex2 not in exonDic2:
				continue
			origExonEntry1=exonDic1[ex1] # org1 
			origExonEntry2=exonDic2[ex2] # org2
			if ex1 not in liftoverExonMappings:
				liftoverExonMappings[ex1]={}
			if ex2 not in liftoverExonMappings[ex1]:
				# first for allLength, second for coding if partCoding or fullCoding, then minMatch levels
				liftoverExonMappings[ex1][ex2]=[0,0,"0","0"] 
			#
		# if oneToTwo is False. now the liftover length is for org2 and so should the origExonEntry1 be
		else: 
			# discard those that has nothing to do with this gene
			if ex2 not in exonDic1:
				continue
			
			origExonEntry2=exonDic1[ex2] # org1
			origExonEntry1=exonDic2[ex1] # org2

			if ex2 not in liftoverExonMappings:
				liftoverExonMappings[ex2]={}
			if ex1 not in liftoverExonMappings[ex2]:
				# first for allLength, second for coding if partCoding or fullCoding, then minMatch levels
				liftoverExonMappings[ex2][ex1]=[0,0,"0","0"] 
				#print ["********************************NEW",ex2,ex1]
			#
			ex1,ex2=ex2,ex1
		#

		origExonLen1=abs(origExonEntry1.basicInfoDic["end_coord"]-origExonEntry1.basicInfoDic["start_coord"])
		origExonCodingLen1=abs(origExonEntry1.codingExon[0]-origExonEntry1.codingExon[1]) # will be zero for nonCoding and full length for fullCoding
		origExonLen2=abs(origExonEntry2.basicInfoDic["end_coord"]-origExonEntry2.basicInfoDic["start_coord"])
		origExonCodingLen2=abs(origExonEntry2.codingExon[0]-origExonEntry2.codingExon[1])

		overlap=max(0, min(int(e1),int(e2))-max(int(s1),int(s2)))
		# entry is reporting overlap between two partCoding exons
		if suffixStr.endswith("partCoding"):
			codingOverlapScore=compute_overlap_score_for_exons(origExonCodingLen1,origExonCodingLen2,liftedOverLenForEx,overlap)
		#	print [oneToTwo, "partCoding", ex1, ex2, origExonCodingLen1,origExonCodingLen2,liftedOverLenForEx,overlap, codingOverlapScore]
		#	print [oneToTwo, "partCoding2", ex1, ex2, s1,e1,origExonEntry2.codingExon]
		else: # entry is reporting overlap between two full length exons (either fullCoding or nonCoding or partCoding)
			allOverlapScore=compute_overlap_score_for_exons(origExonLen1,origExonLen2,liftedOverLenForEx,overlap)
		#	print [oneToTwo,"fullLength", ex1, ex2, origExonLen1,origExonLen2,liftedOverLenForEx,overlap, allOverlapScore]

		# make sure to set coding overlap to the same number as overall overlap for full coding pairs
		#if origExonEntry1.exon_type=="fullCoding" and origExonEntry2.exon_type=="fullCoding":
		#	codingOverlapScore=allOverlapScore

		if liftoverExonMappings[ex1][ex2][1]<codingOverlapScore or \
			(liftoverExonMappings[ex1][ex2][1]==codingOverlapScore and float(liftoverExonMappings[ex1][ex2][3])<float(minMatchStr)):
			liftoverExonMappings[ex1][ex2][1]=codingOverlapScore
			liftoverExonMappings[ex1][ex2][3]=minMatchStr
		if liftoverExonMappings[ex1][ex2][0]<allOverlapScore or \
			(liftoverExonMappings[ex1][ex2][0]==allOverlapScore and float(liftoverExonMappings[ex1][ex2][2])<float(minMatchStr)):
			liftoverExonMappings[ex1][ex2][0]=allOverlapScore
			liftoverExonMappings[ex1][ex2][2]=minMatchStr
	#
	infile.close()
	return 


def parse_liftOver_finalExonMappings_perGenePair(genePairId):
	"""
	Parse the three different types of files for each exon for 
	this gene pair and return the summaries for the next step.
	"""
	global liftoverExonMappings,unmappedExons,nonintersectingExons
	liftoverExonMappings={} 
	unmappedExons={}
	nonintersectingExons={}
	#sys.stderr.write("Parsing liftover final exon mappings from single exon files for both organisms\n\n")

	c=1
	for exonDic in [exonDic1,exonDic2]:
		isOneToTwo = True if c==1 else False #ifelse one liner
		counts={"000": 0, "001": 0, "010": 0, "011": 0, "100": 0, "101": 0, "110": 0, "111": 0}
		#
		for exonId in exonDic:
			infilenameMapped=perExonLiftoverDir+"/org"+str(c)+"/"+exonId+"_mapped.txt"
			infilenameUnmapped=perExonLiftoverDir+"/org"+str(c)+"/"+exonId+"_unmapped.txt"
			infilenameNonintersecting=perExonLiftoverDir+"/org"+str(c)+"/"+exonId+"_nonintersecting.txt"
			
			isM,isU,isN=os.path.isfile(infilenameMapped),os.path.isfile(infilenameUnmapped),os.path.isfile(infilenameNonintersecting)
			#	counts["001"]=counts["001"]+1 # only nonintersecting - ~6k 
			#	counts["010"]=counts["010"]+1 # only unmapped - ~69k
			#	counts["011"]=counts["011"]+1 # unmapped and nonintersecting - ~13k
			#	counts["100"]=counts["100"]+1  # only mapped - ~258k
			#	counts["101"]=counts["101"]+1 # mapped and nonintersecting - ~60k
			#	counts["110"]=counts["110"]+1 # mapped and unmapped - ~83k
			#	counts["111"]=counts["111"]+1 # all three files - ~11k
			
			m,u,n=0,0,0
			if isM:
				m=1
				parse_liftOver_mapped_file_perExon(infilenameMapped,isOneToTwo) # True means it is from org1 to org2
			if isU:
				u=1
				unmappedExons[exonId]={}
				parse_liftOver_unmapped_file_perExon(infilenameUnmapped)
			if isN:
				n=1
				parse_liftOver_nonintersecting_file_perExon(infilenameNonintersecting,isOneToTwo) # True means it is from org1 to org2
			#
			munStr=str(m)+str(u)+str(n)
			counts[munStr]=counts[munStr]+1
		#
		
		print "SummaryPerGene\torg"+str(c),
		for k in ["000", "001", "010", "011", "100", "101", "110", "111"]:
			print ("\t%d" % counts[k]),
		print

		c+=1
	#

	return # return from parse_liftOver_finalExonMappings_perGenePair

def sort_exonDics_byCoordinates(origExonCoords):
	"""
	Order the exons wrt their start coordinates and delete the ones 
	that are duplicates in terms of coordinates. Also report duplicates.
	"""
	duplicateExonsDic={}
	startAndEndCoords={}
	keepDupTrack={}
	toDeleteList=[]

	for exId in origExonCoords:
		dummy,ch,s,e,strand=origExonCoords[exId]
		key=ch+":"+str(s)+"-"+str(e)
		if key not in keepDupTrack:
			keepDupTrack[key]=[]
			startAndEndCoords[exId]=[s,e]
		else:
			toDeleteList.append(exId)
		keepDupTrack[key].append(exId)
	#

	for key in keepDupTrack:
		allExons=keepDupTrack[key]
		if len(allExons)==1:
			continue
		for e in allExons:
			duplicateExonsDic[e]=allExons
	#
	#print keepDupTrack

	strandLast=strand # will use this to make sure all exons sorting is in right order
	if strandLast=="+":
		sortedExons=sorted(startAndEndCoords, key=startAndEndCoords.get) # sort according to start sites then end sites
	else:
		sortedExons=sorted(startAndEndCoords, key=startAndEndCoords.get, reverse=True) # sort according to start sites then end sites
	# 

	# overwrite the first field in the list with exon order, it was redundant anyway. 
	c=0
	for exId in sortedExons:
		origExonCoords[exId][0]=c 
		c+=1
	#
	# delete exons that were not included in sortedExons
	for exId in toDeleteList:
		del origExonCoords[exId]
	#	
	return duplicateExonsDic


def findMatchingsPerGenePair(genePairId):
	global duplicateExons1,duplicateExons2

	# parse exons in the desired format for below functions per this gene pair
	origExonCoords1={}; 
	for exonId in exonDic1:
		ex=exonDic1[exonId]
		ch,s,e,strand="chr"+ex.basicInfoDic["chromosome"],ex.basicInfoDic["start_coord"],ex.basicInfoDic["end_coord"],ex.basicInfoDic["strand"]
		origExonCoords1[exonId]=[exonId,ch,s,e,strand]
	#
	duplicateExons1=sort_exonDics_byCoordinates(origExonCoords1)
	
	origExonCoords2={}
	for exonId in exonDic2:
		ex=exonDic2[exonId]
		ch,s,e,strand="chr"+ex.basicInfoDic["chromosome"],ex.basicInfoDic["start_coord"],ex.basicInfoDic["end_coord"],ex.basicInfoDic["strand"]
		origExonCoords2[exonId]=[exonId,ch,s,e,strand]
	#
	duplicateExons2=sort_exonDics_byCoordinates(origExonCoords2)

	########################################################################################### 
	### NOW origExonCoords1,origExonCoords2 have exons after duplicate removal and in order ###
	### exonDic1,exonDic2 have exons before duplicate removal in no specific order 		###
	########################################################################################### 

	print "*****************************************************************"
	print "Number of exons before and after duplicate removal"
	print "Org1\t%d\t%d\n" % (len(exonDic1),len(origExonCoords1)),
	print "Org2\t%d\t%d\n" % (len(exonDic2),len(origExonCoords2)),
	print "*****************************************************************"

	# After this function call I will have all the information about mapped, unmapped, nonintersecting exons
	# in these three dictionaries: liftoverExonMappings, unmappedExons, nonintersectingExons
	# NOTE THAT these info will be at the exon level before duplicate removal!!! 
	parse_liftOver_finalExonMappings_perGenePair(genePairId)

	# Now decide what class each exon will be assigned to among Mapped, Unmapped, Nonintersecing or left unassigned (not in exonClasses)
	exonClasses1,exonClasses2=classify_exons(genePairId,MAPPED_EXON_THRES)

	# pass this outfile which will be closed in the below function FIXME: uncomment this at some point
	output_transcript_level_matches(genePairId,exonClasses1,exonClasses2)

	return #exonLevelMatchScoresWithDups,transcriptLevelMatches

def classify_exons(genePairId,mappedExonThreshold):
	"""
	Given a threshold on the similarity between two exons to be
	deemed "mapped" (mappedExonThreshold), this function determines
	the exon classes using all the info from parsed liftover files.
	"""
	# open the output file to which the results of the mapping will be written
	outfilename=outdir+"/exonClasses-"+str(MAPPED_EXON_THRES)+".txt"
	outfile=open(outfilename,'w')
	#g1,g2=genePairId.split("-")
	outfile.write("GeneInfo\t%s\n" % geneDic1[g1].get_summary_string())
	outfile.write("GeneInfo\t%s\n" % geneDic2[g2].get_summary_string())
	
	exonClasses1={}; exonClasses2={} 
	## NOW print all that you learned from parsing these three files per exon
	for ex1 in liftoverExonMappings:
		origExonEntry1=exonDic1[ex1]
		ch1,s1="chr"+origExonEntry1.basicInfoDic["chromosome"],origExonEntry1.basicInfoDic["start_coord"]
		e1,strand1=origExonEntry1.basicInfoDic["end_coord"],origExonEntry1.basicInfoDic["strand"]
		e1type=origExonEntry1.exon_type
		for ex2 in liftoverExonMappings[ex1]:
			origExonEntry2=exonDic2[ex2]
			ch2,s2="chr"+origExonEntry2.basicInfoDic["chromosome"],origExonEntry2.basicInfoDic["start_coord"]
			e2,strand2=origExonEntry2.basicInfoDic["end_coord"],origExonEntry2.basicInfoDic["strand"]
			e2type=origExonEntry2.exon_type

			allOverlapScore,codingOverlapScore,minMatchStrAll,minMatchStrCoding=liftoverExonMappings[ex1][ex2]
			if allOverlapScore>=mappedExonThreshold or codingOverlapScore>=mappedExonThreshold:
				if ex1 not in exonClasses1:
					exonClasses1[ex1]=["Type",0,0]
				if ex2 not in exonClasses2:
					exonClasses2[ex2]=["Type",0,0]
				if allOverlapScore>exonClasses1[ex1][1]:
					exonClasses1[ex1][0]="Mapped"; exonClasses1[ex1][1]=allOverlapScore
				if codingOverlapScore>exonClasses1[ex1][2]:
					exonClasses1[ex1][0]="Mapped"; exonClasses1[ex1][2]=codingOverlapScore
				if allOverlapScore>exonClasses2[ex2][1]:
					exonClasses2[ex2][0]="Mapped"; exonClasses2[ex2][1]=allOverlapScore
				if codingOverlapScore>exonClasses2[ex2][2]:
					exonClasses2[ex2][0]="Mapped"; exonClasses2[ex2][2]=codingOverlapScore
				#
			#
			outfile.write("MatchedExonPair\t%s\t%s\t%.3f\t%.3f\t%s\t%s\n" % \
				(origExonEntry1.get_summary_string(),origExonEntry2.get_summary_string(),allOverlapScore,\
				codingOverlapScore,minMatchStrAll,minMatchStrCoding))
			print("%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%.3f\t%.3f\t%s\t%s\tMatchedExonPair\n" \
				% (ch1,s1,e1,strand1,ex1,e1type,ch2,s2,e2,strand2,ex2,e2type,allOverlapScore,codingOverlapScore,minMatchStrAll,minMatchStrCoding)),
		#	
	#

	for ex in nonintersectingExons:
		for chse in nonintersectingExons[ex]:
			minMatchStr,junctionProblem,ch,s,e,strand,sqAcc,sqDon,liftedOverLenForEx=nonintersectingExons[ex][chse]
			if ex in exonDic1:
				origExonEntry=exonDic1[ex]
				# if a mapping not found for this exon with at least mappedExonThreshold and there were no Unmapped entries
				if ex not in exonClasses1:
					exonClasses1[ex]=["Nonintersecing", junctionProblem, minMatchStr]
			else:
				origExonEntry=exonDic2[ex]
				# if a mapping not found for this exon with at least mappedExonThreshold and there were no Unmapped entries
				if ex not in exonClasses2:
					exonClasses2[ex]=["Nonintersecing", junctionProblem, minMatchStr]
			#
			# this doesn't mean that Nonintersecing is the class of this exon
			outfile.write("PotentialSpliceLoss\t%r\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%d\n" % \
				(junctionProblem,minMatchStr,origExonEntry.get_summary_string(),ch,int(s),int(e),strand,sqAcc,sqDon,liftedOverLenForEx))
		#
	#

	for ex in unmappedExons:
		for whyUnmapped in unmappedExons[ex]:
			minMatchStr=unmappedExons[ex][whyUnmapped]
			if ex in exonDic1:
				origExonEntry=exonDic1[ex]
				# if a mapping not found for this exon with at least mappedExonThreshold
				if ex not in exonClasses1:
					exonClasses1[ex]=["Unmapped",whyUnmapped,minMatchStr]
			else:
				origExonEntry=exonDic2[ex]
				# if a mapping not found for this exon with at least mappedExonThreshold
				if ex not in exonClasses2:
					exonClasses2[ex]=["Unmapped",whyUnmapped,minMatchStr]
			#
			# this doesn't mean that UnmappedExon is the class of this exon
			outfile.write("UnmappedExon\t%s\t%s\t%s\n" % (minMatchStr,origExonEntry.get_summary_string(),whyUnmapped))
		#
	#
	print("***************** Exon classes for the first organism **************************\n"),
	for ex in exonDic1:
		origExonEntry=exonDic1[ex]
		cl,why,prm="OTHER",0,0
		if ex in exonClasses1:
			cl,why,prm=exonClasses1[ex][0],exonClasses1[ex][1],exonClasses1[ex][2]
		#
		outfile.write("ExonClass:\t%s\t%r\t%r\t%r\n" % (origExonEntry.get_summary_string(),cl,why,prm))
		print("ExonClass:\t%s\t%r\t%r\t%r\n" % (origExonEntry.get_summary_string(),cl,why,prm)),

	print("***************** Exon classes for the second organism **************************\n"),
	for ex in exonDic2:
		origExonEntry=exonDic2[ex]
		cl,why,prm="OTHER",0,0
		if ex in exonClasses2:
			cl,why,prm=exonClasses2[ex][0],exonClasses2[ex][1],exonClasses2[ex][2]
		#
		outfile.write("ExonClass:\t%s\t%r\t%r\t%r\n" % (origExonEntry.get_summary_string(),cl,why,prm))
		print("ExonClass:\t%s\t%r\t%r\t%r\n" % (origExonEntry.get_summary_string(),cl,why,prm)),
	print("*********************************************************************************\n"),
	outfile.close()
	return exonClasses1,exonClasses2

def dynamicProgramming_per_transcript_pair(t1,t2):

#	def get_score_all(i, j):
#		ex1=elist1[i-1][2]
#		ex2=elist2[j-1][2]
#		allOverlapScore,codingOverlapScore=0,0
#		if ex1 in liftoverExonMappings and ex2 in liftoverExonMappings[ex1]:
#			allOverlapScore,codingOverlapScore,minMatchStrAll,minMatchStrCoding=liftoverExonMappings[ex1][ex2]
#			#
#		#
#		return allOverlapScore

	def get_score_coding(i, j):
		ex1=elist1[i-1][2]
		ex2=elist2[j-1][2]
		e1type=exonDic1[ex1].exon_type
		e2type=exonDic2[ex2].exon_type
		allOverlapScore,codingOverlapScore=0,0
		if ex1 in liftoverExonMappings and ex2 in liftoverExonMappings[ex1]:
			allOverlapScore,codingOverlapScore,minMatchStrAll,minMatchStrCoding=liftoverExonMappings[ex1][ex2]
			if e1type=="fullCoding" and e2type=="fullCoding":
				codingOverlapScore=allOverlapScore
			#
		#
		return codingOverlapScore

	def get_score_all(i, j):
		ex1=elist1[i-1][2]
		ex2=elist2[j-1][2]
		e1type=exonDic1[ex1].exon_type
		e2type=exonDic2[ex2].exon_type
		allOverlapScore,codingOverlapScore=0,0
		if ex1 in liftoverExonMappings and ex2 in liftoverExonMappings[ex1]:
			allOverlapScore,codingOverlapScore,minMatchStrAll,minMatchStrCoding=liftoverExonMappings[ex1][ex2]
			if e1type=="fullCoding" and e2type=="fullCoding":
				codingOverlapScore=allOverlapScore
			#
		#
		return max(codingOverlapScore,allOverlapScore)
	

	elist1=transcriptDic1[t1].exons
	elist2=transcriptDic2[t2].exons
	ecodinglist1=[]
	for ex1 in [x[2] for x in elist1]:
		if exonDic1[ex1].exon_type!="nonCoding": ecodinglist1.append(ex1)
	ecodinglist2=[]
	for ex2 in [x[2] for x in elist2]:
		if exonDic2[ex2].exon_type!="nonCoding": ecodinglist2.append(ex2)
			
	noRows=len(elist1)
	noCols=len(elist2)
	dpTableAll=np.zeros((noRows+1,noCols+1))
	dpTableCoding=np.zeros((noRows+1,noCols+1))
	dppathAll=[]
	dppathCoding=[]

	indel=0 # indel or gap penalty

	for i in xrange(noRows+1):
		dpTableAll[i,0] = i*indel
		dpTableCoding[i,0] = i*indel
	for j in xrange(noCols+1):
		dpTableAll[0,j] = j*indel
		dpTableCoding[0,j] = j*indel

	#noOfAlignedExonPairs=0
	for i in xrange(1, noRows+1):
		for j in xrange(1, noCols+1):
			# first all: which is the max between whole and coding
			dpTableAll[i,j]=max(dpTableAll[i-1, j-1]+ get_score_all(i, j), dpTableAll[i-1, j] + indel, dpTableAll[i, j-1] + indel)
			#if get_score_all(i, j)>0 and dpTableAll[i,j]-dpTableAll[i-1, j-1]==get_score_all(i, j):
			#	dppathAll.append([i-1,j-1,get_score_all(i, j),get_score_coding(i, j)])
			# then coding
			dpTableCoding[i,j]=max(dpTableCoding[i-1, j-1]+ get_score_coding(i, j), dpTableCoding[i-1, j] + indel, dpTableCoding[i, j-1] + indel)
			#if get_score_coding(i, j)>0 and dpTableCoding[i,j]-dpTableCoding[i-1, j-1]==get_score_coding(i, j):
			#	dppathCoding.append([i-1,j-1,get_score_all(i, j),get_score_coding(i, j)])
			
	#
	#print ["DP\t",t1,t2]
	#print dpTableAll
	dpscoreAll=dpTableAll[noRows,noCols]
	print "DPall\t%s\t%s\t%d\t%d\t%.3f\t%.3f\t%s\t%s\t%s\t%s\n" % (t1,t2,len(elist1),len(elist2),dpscoreAll,dpscoreAll/max(len(elist1),len(elist2)),transcriptDic1[t1].basicInfoDic["transcript_name"],transcriptDic2[t2].basicInfoDic["transcript_name"],transcriptDic1[t1].basicInfoDic["transcript_biotype"],transcriptDic2[t2].basicInfoDic["transcript_biotype"]),
	
	dpscoreCoding=dpTableCoding[noRows,noCols]	
	print "DPcoding\t%s\t%s\t%d\t%d\t%.3f\t%.3f\t%s\t%s\t%s\t%s\n" % (t1,t2,len(ecodinglist1),len(ecodinglist2),dpscoreCoding,dpscoreCoding/max(len(ecodinglist1),len(ecodinglist2)),transcriptDic1[t1].basicInfoDic["transcript_name"],transcriptDic2[t2].basicInfoDic["transcript_name"],transcriptDic1[t1].basicInfoDic["transcript_biotype"],transcriptDic2[t2].basicInfoDic["transcript_biotype"]),

	i=noRows; j=noCols
	while i>0 and j>0:
		if get_score_all(i, j)>0 and dpTableAll[i,j]-dpTableAll[i-1, j-1]==get_score_all(i, j):
			dppathAll.append([i-1,j-1,get_score_all(i, j),get_score_coding(i, j)])
			i -= 1; j -= 1;
		elif dpTableAll[i,j]-dpTableAll[i-1, j]==indel:
			i -= 1
		else:
			j -= 1
	#
	dppathAll.reverse()

	i=noRows; j=noCols
	while i>0 and j>0:
		if get_score_coding(i, j)>0 and dpTableCoding[i,j]-dpTableCoding[i-1, j-1]==get_score_coding(i, j):
			dppathCoding.append([i-1,j-1,get_score_all(i, j),get_score_coding(i, j)])
			i -= 1; j -= 1;
		elif dpTableCoding[i,j]-dpTableCoding[i-1, j]==indel:
			i -= 1
		else:
			j -= 1
	#
	dppathCoding.reverse()

	## taking max
	#dpscoreAll=dpscoreAll/max(len(elist1),len(elist2))
	#dpscoreCoding=dpscoreCoding/max(1,len(ecodinglist1),len(ecodinglist2))
	## taking sum
	dpscoreAll=(2*dpscoreAll)/max(1,len(elist1)+len(elist2))
	dpscoreCoding=(2*dpscoreCoding)/max(1,len(ecodinglist1)+len(ecodinglist2))

	return (dpscoreAll,dppathAll,dpscoreCoding,dppathCoding)


def allExonPairs_per_transcript_pair(t1,t2):

	elist1=transcriptDic1[t1].exons
	elist2=transcriptDic2[t2].exons
	ecodinglist1=[]
	for ex1 in [x[2] for x in elist1]:
		if exonDic1[ex1].exon_type!="nonCoding": ecodinglist1.append(ex1)
	ecodinglist2=[]
	for ex2 in [x[2] for x in elist2]:
		if exonDic2[ex2].exon_type!="nonCoding": ecodinglist2.append(ex2)
	
		
	tpairTotals=[0,0,0] # numberTried, numberMatched, totalScore
	tCodingPairTotals=[0,0,0] # numberTried, numberMatched, totalScore
	
	noRows=len(elist1)
	noCols=len(elist2)
	expathAll=[]
	expathCoding=[]


	ei=0
	for ex1 in [x[2] for x in elist1]:
		origExonEntry1=exonDic1[ex1]
		ej=0
		e1type=origExonEntry1.exon_type
		
		for ex2 in [x[2] for x in elist2]:
			origExonEntry2=exonDic2[ex2]
			e2type=origExonEntry2.exon_type
			allOverlapScore,codingOverlapScore=0,0
			if ex1 in liftoverExonMappings and ex2 in liftoverExonMappings[ex1]:
				allOverlapScore,codingOverlapScore,minMatchStrAll,minMatchStrCoding=liftoverExonMappings[ex1][ex2]
			if e1type=="fullCoding" and e2type=="fullCoding":
				codingOverlapScore=allOverlapScore

			tpairTotals[0]+=1
			if max(allOverlapScore,codingOverlapScore)>=MAPPED_EXON_THRES:
				tpairTotals[2]+=max(allOverlapScore,codingOverlapScore)
				expathAll.append([ei,ej,allOverlapScore,codingOverlapScore])
				tpairTotals[1]+=1

			if e1type!="nonCoding" and e2type!="nonCoding":
				tCodingPairTotals[0]+=1
				if codingOverlapScore>=MAPPED_EXON_THRES:
					tCodingPairTotals[2]+=codingOverlapScore
					expathCoding.append([ei,ej,allOverlapScore,codingOverlapScore])
					tCodingPairTotals[1]+=1
			#
			#if allOverlapScore>=MAPPED_EXON_THRES or codingOverlapScore>=MAPPED_EXON_THRES:
			#	print("PerTranscriptExonPair\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%.3f\t%.3f\n" % \
			#		(t1,t2,len(elist1),len(elist2),ei+1,ej+1,ex1,ex2,e1type,e2type,allOverlapScore,codingOverlapScore)),
			#	outfile.write("PerTranscriptExonPair\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%.3f\t%.3f\n" % \
			#		(t1,t2,len(elist1),len(elist2),ei+1,ej+1,ex1,ex2,e1type,e2type,allOverlapScore,codingOverlapScore))
			#
			ej+=1
			#
		ei+=1
	#

	## taking max
	#exscoreAll=(tpairTotals[2])/max(len(elist1),len(elist2))
	#exscoreCoding=(tCodingPairTotals[2])/max(1,max(len(ecodinglist1),len(ecodinglist2)))
	## taking sum
	exscoreAll=(2*tpairTotals[2])/max(1,len(elist1)+len(elist2))
	exscoreCoding=(2*tCodingPairTotals[2])/max(1,len(ecodinglist1)+len(ecodinglist2))

	return (exscoreAll,expathAll,exscoreCoding,expathCoding)


def greedyExonPairs_per_transcript_pair(t1,t2):

	elist1=transcriptDic1[t1].exons
	elist2=transcriptDic2[t2].exons
	ecodinglist1=[]
	for ex1 in [x[2] for x in elist1]:
		if exonDic1[ex1].exon_type!="nonCoding": ecodinglist1.append(ex1)
	ecodinglist2=[]
	for ex2 in [x[2] for x in elist2]:
		if exonDic2[ex2].exon_type!="nonCoding": ecodinglist2.append(ex2)
	
	noRows=len(elist1)
	noCols=len(elist2)

	# create the coding and all similarity matrices for exons
	esimMatrixAll=np.zeros((noRows,noCols))
	esimMatrixCoding=np.zeros((noRows,noCols))

	ei=0
	for ex1 in [x[2] for x in elist1]:
		origExonEntry1=exonDic1[ex1]
		ej=0
		e1type=origExonEntry1.exon_type
		for ex2 in [x[2] for x in elist2]:
			origExonEntry2=exonDic2[ex2]
			e2type=origExonEntry2.exon_type

			allOverlapScore,codingOverlapScore=0,0
			if ex1 in liftoverExonMappings and ex2 in liftoverExonMappings[ex1]:
				allOverlapScore,codingOverlapScore,minMatchStrAll,minMatchStrCoding=liftoverExonMappings[ex1][ex2]
			if e1type=="fullCoding" and e2type=="fullCoding":
				codingOverlapScore=allOverlapScore
	
			esimMatrixAll[ei,ej]=allOverlapScore
			esimMatrixCoding[ei,ej]=codingOverlapScore

			ej+=1
			#
		ei+=1
	#

#	print
#	print esimMatrixAll
#	print esimMatrixCoding
#	print
#	print [x[2] for x in elist1]
#	print [x[2] for x in elist2]
#	print ecodinglist1
#	print ecodinglist2
#	print

	oneToOneExs=[]
	erevDic1={}
	erevDic2={}
	c=0
	for ex1 in [x[2] for x in elist1]:
		erevDic1[c]=ex1
		c+=1 
	c=0
	for ex2 in [x[2] for x in elist2]:
		erevDic2[c]=ex2
		c+=1 

	newEsimMatrixAll=esimMatrixAll.copy()
	newEsimMatrixCoding=esimMatrixCoding.copy()

	newMax=0.001
	conditionCounter=[0,0,0,0,0,0,0,0]
	while sum(sum(newEsimMatrixAll))+sum(sum(newEsimMatrixCoding))>0 and newMax>0 and len(oneToOneExs)<min(len(elist1),len(elist2)):

		# select the max among the two score matrices
		newMaxAll=newEsimMatrixAll.max()
		newMaxCoding=newEsimMatrixCoding.max()
		newMax=max(newMaxAll,newMaxCoding)

		#codingIsMax=True ### FIXME: Implement this later to make sure all vs coding being the max is both handled

		#np.set_printoptions(precision=3)
		#print "**HOW all"
		#print (newEsimMatrixAll)
		#print "**HOW coding"
		#print(newEsimMatrixCoding)

		# choose where the max of the two came from
		if newMaxAll>newMaxCoding:
			xs,ys=np.where(abs(newEsimMatrixAll-newMax)<=0.00000000001)
		else:
			xs,ys=np.where(abs(newEsimMatrixCoding-newMax)<=0.00000000001)

		if len(xs)==0 or len(ys)==0: # Condition 0: something wrong, shouldn't have entered the while loop
			conditionCounter[0]+=1
		if len(xs)==1: # Condition 1: one unique maximum for the maximum similarity (either coding or all)
			x=xs[0]; y=ys[0] 
			oneToOneExs.append([x,y,1,newEsimMatrixCoding[x,y],newEsimMatrixAll[x,y],-1,-1])
			newEsimMatrixCoding[x,:]=0; newEsimMatrixCoding[:,y]=0
			newEsimMatrixAll[x,:]=0; newEsimMatrixAll[:,y]=0
			conditionCounter[1]+=1
		else: # multiple equal maximum similarity (either coding or all)
			allMax=0
			for i in range(len(xs)):
				x=xs[i]; y=ys[i]
				if newEsimMatrixAll[x,y]>allMax:
					allMax=newEsimMatrixAll[x,y]
				#
			#
			axs,ays=np.where(abs(newEsimMatrixAll-allMax)<=0.00000000001)
			if len(axs)==1: # Condition 2: one unique maximum in overall matrix after tie in max (either coding or all)
				x=axs[0]; y=ays[0] 
				oneToOneExs.append([x,y,2,newEsimMatrixCoding[x,y],newEsimMatrixAll[x,y],-1,-1])
				newEsimMatrixCoding[x,:]=0; newEsimMatrixCoding[:,y]=0
				newEsimMatrixAll[x,:]=0; newEsimMatrixAll[:,y]=0
				conditionCounter[2]+=1
			else: # multiple equal coding and overall sim- pick the first and report it
				x=axs[0]; y=ays[0] 
				oneToOneExs.append([x,y,3,newEsimMatrixCoding[x,y],newEsimMatrixAll[x,y],-1,-1])
				newEsimMatrixCoding[x,:]=0; newEsimMatrixCoding[:,y]=0
				newEsimMatrixAll[x,:]=0; newEsimMatrixAll[:,y]=0
				conditionCounter[3]+=1
		#
	#
	#print oneToOneExs

	exscoreAll=0
	exscoreCoding=0
	expathAll=[]
	expathCoding=[]
	for l in oneToOneExs:
		ei,ej,condition,codingScore,allScore,codingExonDiff,allExonDiff=l
		ex1=erevDic1[ei]
		ex2=erevDic2[ej]
		exscoreAll+=allScore
		exscoreCoding+=codingScore
		expathAll.append([ei,ej,allScore,codingScore])
		#
		origExonEntry1=exonDic1[ex1]
		origExonEntry2=exonDic2[ex2]
		e1type=origExonEntry1.exon_type
		e2type=origExonEntry2.exon_type
		if e1type!="nonCoding" and e2type!="nonCoding":
			expathCoding.append([ei,ej,allScore,codingScore])
	#	

	## taking max
	#exscoreAll=(exscoreAll)/max(len(elist1),len(elist2))
	#exscoreCoding=(exscoreCoding)/max(1,max(len(ecodinglist1),len(ecodinglist2)))
	## taking sum
	exscoreAll=(2*exscoreAll)/max(1,len(elist1)+len(elist2))
	exscoreCoding=(2*exscoreCoding)/max(1,len(ecodinglist1)+len(ecodinglist2))

	return (exscoreAll,expathAll,exscoreCoding,expathCoding)


def output_transcript_level_matches(genePairId,exonClasses1,exonClasses2):
	
	outfilename=outdir+"/transcriptLevelMatches-"+str(MAPPED_EXON_THRES)+".txt"
	outfile=open(outfilename,'w')

	# transcript level matching
	tdic1=geneDic1[g1].transcripts
	tdic2=geneDic2[g2].transcripts
	
	noRows=len(tdic1)
	noCols=len(tdic2)
	tsimMatrixAll=np.zeros((noRows,noCols))
	tsimMatrixCoding=np.zeros((noRows,noCols))

	ti=0
	for t1 in tdic1:
		elist1=transcriptDic1[t1].exons
		ecodinglist1=[]
		for ex1 in [x[2] for x in elist1]:
			if exonDic1[ex1].exon_type!="nonCoding": ecodinglist1.append(ex1)

		tj=0
		for t2 in tdic2:
			elist2=transcriptDic2[t2].exons

			scoreAll,pathAll,scoreCoding,pathCoding=greedyExonPairs_per_transcript_pair(t1,t2)
			#scoreAll,pathAll,scoreCoding,pathCoding=dynamicProgramming_per_transcript_pair(t1,t2)
			#scoreAll,pathAll,scoreCoding,pathCoding=allExonPairs_per_transcript_pair(t1,t2)

			tsimMatrixAll[ti,tj]=scoreAll
			tsimMatrixCoding[ti,tj]=scoreCoding

			ecodinglist2=[]
			for ex2 in [x[2] for x in elist2]:
				if exonDic2[ex2].exon_type!="nonCoding": ecodinglist2.append(ex2)

			print ("TranscriptMatches:\t%s\t%s\t%d\t%d\t%d\t%.3f\t%d\t%d\t%d\t%.3f\n" % \
				(transcriptDic1[t1].get_summary_string(),transcriptDic2[t2].get_summary_string(),\
				len(pathAll),len(elist1),len(elist2),tsimMatrixAll[ti,tj],\
				len(pathCoding),len(ecodinglist1),len(ecodinglist2),tsimMatrixCoding[ti,tj])),

			outfile.write("TranscriptMatches:\t%s\t%s\t%d\t%d\t%d\t%.3f\t%d\t%d\t%d\t%.3f\n" % \
				(transcriptDic1[t1].get_summary_string(),transcriptDic2[t2].get_summary_string(),\
				len(pathAll),len(elist1),len(elist2),tsimMatrixAll[ti,tj],\
				len(pathCoding),len(ecodinglist1),len(ecodinglist2),tsimMatrixCoding[ti,tj])),

#			print ("TranscriptMatches:\t%s\t%s\t%d\t%d\t%.3f\t%d\t%d\t%.3f\t%d\t%d\t%.3f\t%d\t%d\t%.3f\n" % \
#				(transcriptDic1[t1].get_summary_string(),transcriptDic2[t2].get_summary_string(),\
#				tpairTotals[0],tpairTotals[1],tpairTotals[2],len(elist1),len(elist2),tsimMatrixAll[ti,tj],\
#				tCodingPairTotals[0],tCodingPairTotals[1],tCodingPairTotals[2],len(ecodinglist1),len(ecodinglist2),tsimMatrixCoding[ti,tj])),
#			outfile.write("TranscriptMatches:\t%s\t%s\t%d\t%d\t%.3f\t%d\t%d\t%.3f\t%d\t%d\t%.3f\t%d\t%d\t%.3f\n" % \
#				(transcriptDic1[t1].get_summary_string(),transcriptDic2[t2].get_summary_string(),\
#				tpairTotals[0],tpairTotals[1],tpairTotals[2],len(elist1),len(elist2),tsimMatrixAll[ti,tj],\
#				tCodingPairTotals[0],tCodingPairTotals[1],tCodingPairTotals[2],len(ecodinglist1),len(ecodinglist2),tsimMatrixCoding[ti,tj]))

			for l in pathAll:
				ei,ej,allOverlapScore,codingOverlapScore=l[0],l[1],l[2],l[3]
				ex1=elist1[ei][2]
				ex2=elist2[ej][2]
				e1type=exonDic1[ex1].exon_type
				e2type=exonDic2[ex2].exon_type
				print("PerTranscriptExonPair\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%.3f\t%.3f\n" % \
					(t1,t2,len(elist1),len(elist2),ei+1,ej+1,ex1,ex2,e1type,e2type,allOverlapScore,codingOverlapScore)),
				outfile.write("PerTranscriptExonPair\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%.3f\t%.3f\n" % \
					(t1,t2,len(elist1),len(elist2),ei+1,ej+1,ex1,ex2,e1type,e2type,allOverlapScore,codingOverlapScore))

			tj+=1

		#
		ti+=1
	#return

	trevDic1={}
	trevDic2={}
	c=0
	for t in tdic1:
		trevDic1[c]=t
		c+=1 
	c=0
	for t in tdic2:
		trevDic2[c]=t
		c+=1 

	oneToOneTrs=greedily_map_transcripts_with_tie_breaks(tdic1,tdic2,tsimMatrixAll,tsimMatrixCoding)

	outfilename=perGenePairPickleDir+"/"+genePairId+"/org1-ucsc-"+str(MAPPED_EXON_THRES)+".bed"
	outfile1=open(outfilename,'w')
	gEntry=geneDic1[g1]
	outfile1.write("browser position chr%s:%d-%d\n" % (gEntry.basicInfoDic["chromosome"],gEntry.basicInfoDic["start_coord"],gEntry.basicInfoDic["end_coord"]))
	outfile1.write("browser hide all\n")
	outfile1.write("track name=\"%s\" description=\"%s-%s gene pair\" visibility=2 color=0,128,0 useScore=1\n" % (gEntry.basicInfoDic["gene_name"],g1,g2))

	outfilename=perGenePairPickleDir+"/"+genePairId+"/org2-ucsc-"+str(MAPPED_EXON_THRES)+".bed"
	outfile2=open(outfilename,'w')
	gEntry=geneDic2[g2]
	outfile2.write("browser position chr%s:%d-%d\n" % (gEntry.basicInfoDic["chromosome"],gEntry.basicInfoDic["start_coord"],gEntry.basicInfoDic["end_coord"]))
	outfile2.write("browser hide all\n")
	outfile2.write("track name=\"%s\" description=\"%s-%s gene pair\" visibility=2 color=0,128,0 useScore=1\n" % (gEntry.basicInfoDic["gene_name"],g2,g1))


	for l in oneToOneTrs:
		ti,tj,condition,codingScore,allScore,codingExonDiff,allExonDiff=l
		t1Entry=transcriptDic1[trevDic1[ti]]
		t2Entry=transcriptDic2[trevDic2[tj]]
		print ("OneToOneTrMatches:\t%s\t%s\t%d\t%.3f\t%.3f\t%d\t%d\n" % \
			(t1Entry.get_summary_string(),t2Entry.get_summary_string(),condition,codingScore,allScore,codingExonDiff,allExonDiff)),
		outfile.write("OneToOneTrMatches:\t%s\t%s\t%d\t%.3f\t%.3f\t%d\t%d\n" % \
			(t1Entry.get_summary_string(),t2Entry.get_summary_string(),condition,codingScore,allScore,codingExonDiff,allExonDiff))


		############ do this for t1Entry ########################
		ch,st,en="chr"+t1Entry.basicInfoDic["chromosome"],t1Entry.basicInfoDic["start_coord"],t1Entry.basicInfoDic["end_coord"]
		name=t1Entry.basicInfoDic["transcript_name"]
		color=900
		strand=t1Entry.basicInfoDic["strand"]
		thickSt,thickEn=ch,st # for now just make everything coding
		rgbset=0 # (r,g,b) triplet if wanted to use
		noOfExons=len(t1Entry.exons)
		blockSizes=""
		blockStarts=""
		if strand=="+":
			for exEntry in t1Entry.exons:
				eSt,eEn=exEntry[0],exEntry[1]
				blockSizes=blockSizes+","+str(abs(eEn-eSt))
				blockStarts=blockStarts+","+str(eSt-st) # exon offsets
			#
			blockSizes="0"+blockSizes+",1" 
			blockStarts="0"+blockStarts+","+str(abs(en-st)-1)
		else:
			for exEntry in t1Entry.exons:
				eSt,eEn=exEntry[0],exEntry[1]
				blockSizes=str(abs(eEn-eSt))+","+blockSizes
				blockStarts=str(eSt-st)+","+blockStarts # exon offsets
			#
			blockSizes="0,"+blockSizes+"1" 
			blockStarts="0,"+blockStarts+str(abs(en-st)-1) 
		#

		minSt,maxEn=en,st # initialize so that min, max works fine
		for exEntry in t1Entry.codingExons:
			eSt,eEn=exEntry[0],exEntry[1]
			if eSt<minSt and eSt>-1:
				minSt=eSt
			if eEn>maxEn and eEn>-1:
				maxEn=eEn
		#
		thickSt,thickEn=minSt,maxEn
		# write the original transcript for org1
		outfile1.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%s\n" % \
			(ch,st,en,name,color,strand,thickSt,thickEn,rgbset,noOfExons+2,blockSizes,blockStarts))
		############ do this for t1Entry - DONE ########################

		########### now find the liftover partners of each exon to print into outfile1 #########################
		liftedOverExons={}
		for ex1Entry in t1Entry.exons:
			exId1=ex1Entry[2]
			if exId1 in exonClasses1 and exonClasses1[exId1][0]=="Mapped":
				liftedOverExons[exId1]=1 # has a liftedover version - conserved
			#
		#
		name=t2Entry.basicInfoDic["transcript_name"] # get the name from the other org
		color=300
		noOfExons=0 
		blockSizes=""
		blockStarts=""
		if strand=="+":
			for exEntry in t1Entry.exons:
				eSt,eEn,eId=exEntry
				if eId in liftedOverExons:
					blockSizes=blockSizes+","+str(abs(eEn-eSt))
					blockStarts=blockStarts+","+str(eSt-st) # exon offsets
					noOfExons+=1
			#
			blockSizes="0"+blockSizes+",1" 
			blockStarts="0"+blockStarts+","+str(abs(en-st)-1)
		else:
			for exEntry in t1Entry.exons:
				eSt,eEn,eId=exEntry
				if eId in liftedOverExons:
					blockSizes=str(abs(eEn-eSt))+","+blockSizes
					blockStarts=str(eSt-st)+","+blockStarts # exon offsets
					noOfExons+=1
			#
			blockSizes="0,"+blockSizes+"1" 
			blockStarts="0,"+blockStarts+str(abs(en-st)-1) 
		#
		# write the original transcript for org1's correspondence in org2
		outfile1.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%s\n" % \
			(ch,st,en,name,color,strand,thickSt,thickEn,rgbset,noOfExons+2,blockSizes,blockStarts))
		######################################################################################################

		############ do this for t2Entry ########################
		ch,st,en="chr"+t2Entry.basicInfoDic["chromosome"],t2Entry.basicInfoDic["start_coord"],t2Entry.basicInfoDic["end_coord"]
		name=t2Entry.basicInfoDic["transcript_name"]
		color=900
		strand=t2Entry.basicInfoDic["strand"]
		thickSt,thickEn=ch,st # for now just make everything coding
		rgbset=0 # (r,g,b) triplet if wanted to use
		noOfExons=len(t2Entry.exons)
		blockSizes=""
		blockStarts=""
		if strand=="+":
			for exEntry in t2Entry.exons:
				eSt,eEn=exEntry[0],exEntry[1]
				blockSizes=blockSizes+","+str(abs(eEn-eSt))
				blockStarts=blockStarts+","+str(eSt-st) # exon offsets
			#
			blockSizes="0"+blockSizes+",1" 
			blockStarts="0"+blockStarts+","+str(abs(en-st)-1)
		else:
			for exEntry in t2Entry.exons:
				eSt,eEn=exEntry[0],exEntry[1]
				blockSizes=str(abs(eEn-eSt))+","+blockSizes
				blockStarts=str(eSt-st)+","+blockStarts # exon offsets
			#
			blockSizes="0,"+blockSizes+"1" 
			blockStarts="0,"+blockStarts+str(abs(en-st)-1) 
		#

		minSt,maxEn=en,st # initialize so that min, max works fine
		for exEntry in t2Entry.codingExons:
			eSt,eEn=exEntry[0],exEntry[1]
			if eSt<minSt and eSt>-1:
				minSt=eSt
			if eEn>maxEn and eEn>-1:
				maxEn=eEn
		#
		thickSt,thickEn=minSt,maxEn
		# write the original transcript for org2
		outfile2.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%s\n" % \
			(ch,st,en,name,color,strand,thickSt,thickEn,rgbset,noOfExons+2,blockSizes,blockStarts))
		############ do this for t2Entry - DONE ########################

		########### now find the liftover partners of each exon to print into outfile2 #########################
		liftedOverExons={}
		for ex2Entry in t2Entry.exons:
			exId2=ex2Entry[2]
			if exId2 in exonClasses2 and exonClasses2[exId2][0]=="Mapped":
				liftedOverExons[exId2]=1 # has a liftedover version - conserved
			#
		#
		name=t1Entry.basicInfoDic["transcript_name"] # get the name from the other org
		color=300
		noOfExons=0 
		blockSizes=""
		blockStarts=""
		if strand=="+":
			for exEntry in t2Entry.exons:
				eSt,eEn,eId=exEntry
				if eId in liftedOverExons:
					blockSizes=blockSizes+","+str(abs(eEn-eSt))
					blockStarts=blockStarts+","+str(eSt-st) # exon offsets
					noOfExons+=1
			#
			blockSizes="0"+blockSizes+",1" 
			blockStarts="0"+blockStarts+","+str(abs(en-st)-1)
		else:
			for exEntry in t2Entry.exons:
				eSt,eEn,eId=exEntry
				if eId in liftedOverExons:
					blockSizes=str(abs(eEn-eSt))+","+blockSizes
					blockStarts=str(eSt-st)+","+blockStarts # exon offsets
					noOfExons+=1
			#
			blockSizes="0,"+blockSizes+"1" 
			blockStarts="0,"+blockStarts+str(abs(en-st)-1) 
		#
		# write the original transcript for org2's correspondence in org1
		outfile2.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%s\n" % \
			(ch,st,en,name,color,strand,thickSt,thickEn,rgbset,noOfExons+2,blockSizes,blockStarts))
		######################################################################################################
		#chr22 2000 7000 itemB 200 - 2200 6950 0 4 433,100,550,1500 0,500,2000,3500
		# NEW STUFF - END
	####

	#print [t for t in tdic1]
	#print [t for t in tdic2]
	#print tsimMatrixAll
	#print tsimMatrixCoding
	outfile.close()
	return

def greedily_map_transcripts_with_tie_breaks(tdic1,tdic2,tsimMatrixAll,tsimMatrixCoding):
	"""
	1. Pick the highest scoring match in the coding matrix 
		- if it is unique then add this matching and set the corresponding column and row to zeros
	2. If there is a tie then look at the all matrix
	3. If there is a tie then look at the difference between the number of coding exons
	4. If there is a tie then look at the difference between the number of exons
	## IMPLEMENT BELOW CONDITIONS AS WELL - FIXME
	5. If there is a tie then look at the sum of the difference between the lengths of coding exons
	6. If still tied then pick one randomly and report the situation
	"""
	
	oneToOneTrs=[]
	trevDic1={}
	trevDic2={}
	c=0
	for t in tdic1:
		trevDic1[c]=t
		c+=1 
	c=0
	for t in tdic2:
		trevDic2[c]=t
		c+=1 

	newTsimMatrixAll=tsimMatrixAll.copy()
	newTsimMatrixCoding=tsimMatrixCoding.copy()

	#print newTsimMatrixAll
	#print newTsimMatrixCoding

	noTr1s=len(tdic1)
	noTr2s=len(tdic2)
	newMax=0.001
	conditionCounter=[0,0,0,0,0,0,0,0]
	#FIXME: I dropped the condition that ensured not too many mappings to allow ties to be reported
	#while sum(sum(newTsimMatrixCoding))>0 and newMax>0 and len(oneToOneTrs)<min(noTr1s,noTr2s):
	while sum(sum(newTsimMatrixAll))+sum(sum(newTsimMatrixCoding))>0 and newMax>0:

		newMaxAll=newTsimMatrixAll.max()
		newMaxCoding=newTsimMatrixCoding.max()
		newMax=max(newMaxCoding,newMaxAll)
		if newMaxAll<0.2: # FIXME: threshold choice but should be better to avoid very small scores
			break
		if newMaxAll>newMaxCoding:
			xs,ys=np.where(abs(newTsimMatrixAll-newMax)<=0.00000000001)
#			print ["All>Coding",newMax,list(xs),list(ys)]
		else:
			xs,ys=np.where(abs(newTsimMatrixCoding-newMax)<=0.00000000001)
#			print ["All<Coding",newMax,list(xs),list(ys)]
		#

		if len(xs)==0 or len(ys)==0: # Condition 0: something wrong, shouldn't have entered the while loop
			conditionCounter[0]+=1
		if len(xs)==1: # Condition 1: one unique maximum 
			x=xs[0]; y=ys[0] 
			oneToOneTrs.append([x,y,1,newTsimMatrixCoding[x,y],newTsimMatrixAll[x,y],-1,-1])
			newTsimMatrixCoding[x,:]=0; newTsimMatrixCoding[:,y]=0
			newTsimMatrixAll[x,:]=0; newTsimMatrixAll[:,y]=0
			conditionCounter[1]+=1
#			print ["conditionCounter: 1",x,y]
### FIXME!!!!!!!!!!
#		elif len(xs)>1: # This automatically makes anything below unreachable
#			x=xs[0]; y=ys[0] 
#			oneToOneTrs.append([x,y,1,newTsimMatrixCoding[x,y],newTsimMatrixAll[x,y],-1,-1])
#			newTsimMatrixCoding[x,y]=0; newTsimMatrixCoding[x,y]=0
#			newTsimMatrixAll[x,y]=0; newTsimMatrixAll[x,y]=0
#			print 
#			print 
#			print newTsimMatrixCoding
#			print 
#			print 
#			conditionCounter[1]+=1
#
		else: # multiple equal coding sim
#			print "conditionCounter: 2+"
			allMax=0
			axs,ays=[],[]
			# find the max first then find th ones among the ones from the prev step that has that max
			for i in range(len(xs)):
				x=xs[i]; y=ys[i]
				if newTsimMatrixAll[x,y]>allMax:
					allMax=newTsimMatrixAll[x,y]
				#
			#
			for i in range(len(xs)):
				x=xs[i]; y=ys[i]
				if abs(newTsimMatrixAll[x,y]-allMax)<=0.00000000001:
					axs.append(x); ays.append(y)
				#
#				print ["conditionCounter: 2+",x,y,newTsimMatrixAll[x,y],newTsimMatrixCoding[x,y],allMax]
			#
			#axs,ays=np.where(abs(newTsimMatrixAll-allMax)<=0.00000000001) # THIS WAS WRONG because it considered all elements not just the prev subset
			if len(axs)==1: # Condition 2: one unique maximum in all matrix after tie in coding
				x=axs[0]; y=ays[0] 
				oneToOneTrs.append([x,y,2,newTsimMatrixCoding[x,y],newTsimMatrixAll[x,y],-1,-1])
				newTsimMatrixCoding[x,:]=0; newTsimMatrixCoding[:,y]=0
				newTsimMatrixAll[x,:]=0; newTsimMatrixAll[:,y]=0
				conditionCounter[2]+=1
#				print ["conditionCounter: 2",x,y]
			else: # multiple equal coding and overall sim
#				print "conditionCounter: 3+"
#				print len(axs)
				minCodingExonDiffList=[]
				for i in range(len(axs)):
					x=axs[i]; y=ays[i]
					t1=trevDic1[x]; t2=trevDic2[y]; 
					elist1=transcriptDic1[t1].exons
					elist2=transcriptDic2[t2].exons
					ex1Counts={"fullCoding":0, "partCoding":0, "nonCoding":0}
					ex2Counts={"fullCoding":0, "partCoding":0, "nonCoding":0}
					for ex1 in [a[2] for a in elist1]:
						ex1Counts[exonDic1[ex1].exon_type]=ex1Counts[exonDic1[ex1].exon_type]+1
					for ex2 in [a[2] for a in elist2]:
						ex2Counts[exonDic2[ex2].exon_type]=ex2Counts[exonDic2[ex2].exon_type]+1
					#
					tempCodingDiff=abs((ex1Counts["fullCoding"]+ex1Counts["partCoding"])-(ex2Counts["fullCoding"]+ex2Counts["partCoding"]))
					minCodingExonDiffList.append(tempCodingDiff)
					#
					#
#					print ["conditionCounter: 3+",x,y,newTsimMatrixAll[x,y],newTsimMatrixCoding[x,y]]
				#

				minCodingExonDiffList=np.array(minCodingExonDiffList)
				minCodingExonDiff=np.min(minCodingExonDiffList)
				ncis=np.where(minCodingExonDiffList==minCodingExonDiff)[0] # one-sided index

				if len(ncis)==1: # Condition 3: one unique minimum difference in the number of coding exons
					x=axs[ncis[0]]; y=ays[ncis[0]] 
					oneToOneTrs.append([x,y,3,newTsimMatrixCoding[x,y],newTsimMatrixAll[x,y],minCodingExonDiff,-1])
					newTsimMatrixCoding[x,:]=0; newTsimMatrixCoding[:,y]=0
					newTsimMatrixAll[x,:]=0; newTsimMatrixAll[:,y]=0
					conditionCounter[3]+=1
#					print ["conditionCounter: 3",x,y]
				else: # multiple equal coding/ and overall sim and overall coding count
#					print "conditionCounter: 4+"
					minAllExonDiffList=[]
					for i in ncis:
						x=axs[i]; y=ays[i]
						t1=trevDic1[x]; t2=trevDic2[y]; 
						elist1=transcriptDic1[t1].exons
						elist2=transcriptDic2[t2].exons
						tempAllDiff=abs(len(elist1)-len(elist2))
						minAllExonDiffList.append(tempAllDiff)
					#
					minAllExonDiffList=np.array(minAllExonDiffList)
					minAllExonDiff=np.min(minAllExonDiffList)
					nais=np.where(minAllExonDiffList==minAllExonDiff)[0] # one-sided index

					#if len(nais)<1:
					#	print "Not sure how"
						
					if len(nais)==1: # Condition 4: one unique minimum difference in the number of all exons
						x=axs[ncis[nais[0]]]; y=ays[ncis[nais[0]]]
#						print [nais[0],ncis[nais[0]],axs[ncis[nais[0]]],ays[ncis[nais[0]]]]
#						print "conditionCounter: 4"
#						print [x,y,newTsimMatrixCoding[x,y],newTsimMatrixAll[x,y]]
#						print newTsimMatrixCoding
#						print newTsimMatrixAll
						oneToOneTrs.append([x,y,4,newTsimMatrixCoding[x,y],newTsimMatrixAll[x,y],minCodingExonDiff,minAllExonDiff])
						newTsimMatrixCoding[x,:]=0; newTsimMatrixCoding[:,y]=0
						newTsimMatrixAll[x,:]=0; newTsimMatrixAll[:,y]=0
						conditionCounter[4]+=1
#						print ["conditionCounter: 4+",x,y,newTsimMatrixAll[x,y],newTsimMatrixCoding[x,y]]
					else: # multiple equal coding/ and overall sim and overall coding count
	#					print "conditionCounter: 5+"
						minCodingLenDiffList=[]
						for i in nais:
							x=axs[ncis[i]]; y=ays[ncis[i]]
							t1=trevDic1[x]; t2=trevDic2[y]; 
							elist1=transcriptDic1[t1].codingExons
							elist2=transcriptDic2[t2].codingExons
							tlen1=0; tlen2=0
							for stC,enC in elist1:
								tlen1+=abs(int(stC)-int(enC))
							for stC,enC in elist2:
								tlen2+=abs(int(stC)-int(enC))
							tempLenDiff=abs(tlen1-tlen2)
							minCodingLenDiffList.append(tempLenDiff)
						#
						minCodingLenDiffList=np.array(minCodingLenDiffList)
						minCodingLenDiff=np.min(minCodingLenDiffList)
						lcis=np.where(minCodingLenDiffList==minCodingLenDiff)[0] # one-sided index

						#if len(lcis)<1:
						#	print "Not sure how"
							
						if len(lcis)==1: # Condition 4: one unique minimum difference in the number of all exons
							x=axs[ncis[nais[lcis[0]]]]; y=ays[ncis[nais[lcis[0]]]]
	#						print "conditionCounter: 5"
							oneToOneTrs.append([x,y,5,newTsimMatrixCoding[x,y],newTsimMatrixAll[x,y],minCodingExonDiff,minAllExonDiff])
							newTsimMatrixCoding[x,:]=0; newTsimMatrixCoding[:,y]=0
							newTsimMatrixAll[x,:]=0; newTsimMatrixAll[:,y]=0
							conditionCounter[5]+=1
						else: # pick and report all of them
	#						print "conditionCounter: 5"
							for i in lcis:
								#print ["How",i,nais[i],]
								x=axs[ncis[nais[i]]]; y=ays[ncis[nais[i]]]
								oneToOneTrs.append([x,y,6,newTsimMatrixCoding[x,y],newTsimMatrixAll[x,y],minCodingExonDiff,minAllExonDiff])
							for i in lcis:
								x=axs[ncis[nais[i]]]; y=ays[ncis[nais[i]]]
								newTsimMatrixCoding[x,:]=0; newTsimMatrixCoding[:,y]=0
								newTsimMatrixAll[x,:]=0; newTsimMatrixAll[:,y]=0
							conditionCounter[6]+=1
					#
			#
		#
	#
	return oneToOneTrs

#   Testing functionalities
def main(argv):
	
	## FIXME: parameter choice
	global MAPPED_EXON_THRES
	MAPPED_EXON_THRES=0.9
	#MAPPED_EXON_THRES=float(argv[1])
	#print MAPPED_EXON_THRES

	global indir,orgId1,orgId2,outdir
	orgId1=argv[1] 
	orgId2=argv[2] 
	indir=ExTraMapperPath+"/preprocess/data/"+orgId1+"-"+orgId2
	
	global ensemblDir,genomedataDir
	global GTFsummaryDir,perGenePairPickleDir,liftOverFilesDir,perExonLiftoverDir

	ensemblDir=indir+"/ensemblDownloads"
	genomedataDir=indir+"/genomedataArchives"
	GTFsummaryDir=indir+"/GTFsummaries"
	perGenePairPickleDir=indir+"/perGenePairPickledInfo" # directory with one picked filed for gene pair information dictionaries
	liftOverFilesDir=indir+"/liftoverRelatedFiles"
	perExonLiftoverDir=indir+"/perExonLiftoverCoords" # directory for one file per exon liftover mappings

	#genePairId="ENSG00000000003-ENSMUSG00000067377"
	#genePairId="ENSG00000197410-ENSMUSG00000102692" # example of a chimeric mouse transcript that maps to two separate human transcripts
	#genePairId="ENSG00000000457-ENSMUSG00000026584" # example that I tested the duplicate removal part on
	#genePairId="ENSG00000000938-ENSMUSG00000028874"
	#genePairId="ENSG00000005020-ENSMUSG00000059182"
	#genePairId="ENSG00000005007-ENSMUSG00000058301"
	#genePairId="ENSG00000003096-ENSMUSG00000036782"
	#genePairId="ENSG00000002745-ENSMUSG00000029671"
	#genePairId="ENSG00000001629-ENSMUSG00000040351"
	#genePairId="ENSG00000002933-ENSMUSG00000023367"
	#genePairId="ENSG00000124203-ENSMUSG00000050600" # one-to-one transcript matching with 5 exons on each side
	#genePairId="ENSG00000146085-ENSMUSG00000023921" # one-to-one transcript matching with 13 exons on each side
	## THIS IS THE EXAMPLE that pointed out the importance of directions in the adjacency graphs ##
	###############################################################################################
	## THIS IS THE EXAMPLE that pointed out that some very small scores included heuristically in the matchings is not good
	#genePairId="ENSG00000258052-ENSMUSG00000064181" # one-to-one transcript matching with 10 exons on one side 5 on the other
	###############################################################################################
	#genePairId="ENSG00000130368-ENSMUSG00000068037" # this is 1 to 20 exon example that used to crash the code since adjacency matrices were all zeros
	#genePairId="ENSG00000121743-ENSMUSG00000048582" # this is a 2 to 2 exon example that used to crash the code since simMatrix was all zeros
	#genePairId="ENSG00000197410-ENSMUSG00000102692" # example of a chimeric mouse transcript that maps to two separate human transcripts
	#genePairId="ENSG00000146085-ENSMUSG00000023921" # one-to-one transcript matching with 13 exons on each side
	#genePairId="ENSG00000185313-ENSMUSG00000034533" # one-to-one transcript matching with 27 exons on one side 28 on the other 
	#genePairId="ENSG00000078900-ENSMUSG00000029026" # TP73 example for Auditi's request
	#genePairId="ENSG00000004948-ENSMUSG00000023964" # CALCR gene with deleted exons from human to mouse
	#genePairId="ENSG00000173769-ENSMUSG00000094985" # TOPAZ1 gene with 20 exons and 1 transcript in each. 2nd exons say deleted but they blast.
	#genePairId="ENSG00000093167-ENSMUSG00000032497" # an exon loss event in mouse due to splice junctions


	#genePairId=argv[1] 
	#findMatchingsPerGenePair(genePairId,perGenePairPickleDir,perExonLiftoverDir)
	#return

	global bigAcceptorDonorSiteDic,bigAcceptorSiteDic,bigDonorSiteDic
	bigAcceptorDonorSiteDic={}
	bigAcceptorSiteDic={}
	bigDonorSiteDic={}
	genePairIds=sorted(os.listdir(perGenePairPickleDir))

	#startIndex=int(argv[1])
	#endIndex=int(argv[2])

	#genePairIds=["ENSG00000140416-ENSMUSG00000032366"]
	#genePairIds=["ENSG00000197410-ENSMUSG00000102692"]
	#genePairIds=["ENSG00000185313-ENSMUSG00000034533"]
	#genePairIds=["ENSG00000173769-ENSMUSG00000094985"]
	#genePairIds=["ENSG00000000457-ENSMUSG00000026584"]
	#genePairIds=["ENSG00000197410-ENSMUSG00000102692"]
	#genePairIds=["ENSG00000146085-ENSMUSG00000023921"]
	#genePairIds=["ENSG00000078900-ENSMUSG00000029026"]
	#genePairIds=["ENSG00000001630-ENSMUSG00000001467"]

	#genePairIds=["ENSG00000140416-ENSMUSG00000032366"] # TPM1 vs Tpm1 
	#genePairIds=["ENSG00000078900-ENSMUSG00000029026"] # TP73
	#genePairIds=["ENSG00000049768-ENSMUSG00000039521"] # FOXP3
#	genePairIds=["ENSG00000157764-ENSMUSG00000002413"] # BRAF
	#genePairIds=["ENSG00000143315-ENSMUSG00000050229"] # 1 exon gene pair
	#genePairIds=["ENSG00000141510-ENSMUSG00000059552"] # TP53

	c=0
	#for genePairId in genePairIds[startIndex:endIndex]:
	#for genePairId in genePairIds[5850:]:
	for genePairId in genePairIds[:10]:
		sys.stderr.write("Finding exon mappings for gene pair number %d\t%s\n" % (c,genePairId))
		global geneDic1,transcriptDic1,exonDic1
		global geneDic2,transcriptDic2,exonDic2
		global g1,g2
		g1,g2=genePairId.split("-")
		# create the output directory
		outdir=ExTraMapperPath+"/output/"+orgId1+"-"+orgId2+"/"+g1+"-"+g2
		if not os.path.exists(outdir): os.makedirs(outdir)

		# reload the pickled dictionaries for this gene pair
		infilename=perGenePairPickleDir+"/"+genePairId+"/org1.pickledDictionaries"
		geneDic1,transcriptDic1,exonDic1=pickle.load(open(infilename,"rb"))
		infilename=perGenePairPickleDir+"/"+genePairId+"/org2.pickledDictionaries"
		geneDic2,transcriptDic2,exonDic2=pickle.load(open(infilename,"rb"))
		origGeneEntry1=geneDic1[g1]; origGeneEntry2=geneDic2[g2]
		# only the ortho pairs with both sides protein_coding

		if origGeneEntry1.basicInfoDic["gene_biotype"]=="protein_coding" and origGeneEntry2.basicInfoDic["gene_biotype"]=="protein_coding":
			findMatchingsPerGenePair(genePairId)
		c+=1
	#
	
	return # from main
	
if __name__ == "__main__":
	main(sys.argv)

