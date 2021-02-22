#!/usr/bin/env python

## ExtraMapper python 3 and later version 
## Original version written by Ferhat Ay
## Converted by Abhijit Chakraborty


import argparse
import sys
import os
import string
import operator
import _pickle as pickle

complement = str.maketrans('atcgn', 'tagcn')

# Reads from exported environment variable
ExTraMapperPath=os.environ['EXTRAMAPPER_DIR']
sys.path.append(ExTraMapperPath+"/scripts")
from ensemblUtils import *

############## Functions ##############

def write_out_bedFiles_for_UCSCbrowser(genePairId,exonClasses1,exonClasses2,oneToOneTrs,trevDic1,trevDic2):

	outfilename=outdir+"/org1-ucsc-"+str(MAPPED_EXON_THRES)+".bed"
	outfile1=open(outfilename,'w')
	gEntry=geneDic1[g1]
	outfile1.write("browser position chr%s:%d-%d\n" % (gEntry.basicInfoDic["chromosome"],gEntry.basicInfoDic["start_coord"],gEntry.basicInfoDic["end_coord"]))
	outfile1.write("browser hide all\n")
	outfile1.write("track name=\"%s\" description=\"%s-%s gene pair\" visibility=2 color=0,128,0 useScore=1\n" % (gEntry.basicInfoDic["gene_name"],g1,g2))
	print ("\nWriting UCSC browser bed output for org1 into file:\n %s" % outfilename)


	outfilename=outdir+"/org2-ucsc-"+str(MAPPED_EXON_THRES)+".bed"
	outfile2=open(outfilename,'w')
	gEntry=geneDic2[g2]
	outfile2.write("browser position chr%s:%d-%d\n" % (gEntry.basicInfoDic["chromosome"],gEntry.basicInfoDic["start_coord"],gEntry.basicInfoDic["end_coord"]))
	outfile2.write("browser hide all\n")
	outfile2.write("track name=\"%s\" description=\"%s-%s gene pair\" visibility=2 color=0,128,0 useScore=1\n" % (gEntry.basicInfoDic["gene_name"],g2,g1))
	print ("Writing UCSC browser bed output for org2 into file:\n %s" % outfilename)

	for l in oneToOneTrs:
		ti,tj,condition,codingScore,allScore,codingExonDiff,allExonDiff=l
		t1Entry=transcriptDic1[trevDic1[ti]]
		t2Entry=transcriptDic2[trevDic2[tj]]
		#print ("OneToOneTrMatches:\t%s\t%s\t%d\t%.3f\t%.3f\t%d\t%d\n" % \
		#	(t1Entry.get_summary_string(),t2Entry.get_summary_string(),condition,codingScore,allScore,codingExonDiff,allExonDiff)),

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
	####

	outfile1.close()
	outfile2.close()
	return

def greedily_map_transcripts_with_tie_breaks(tdic1,tdic2,tsimMatrixAll,tsimMatrixCoding,trevDic1,trevDic2):
	"""
	At each step/condition add the mapping to the overall list and set the corresponding column and row to zeros
	1. Pick the highest scoring match in the coding matrix 
	2. If there is a tie then look at the all matrix
	3. If there is a tie then look at the difference between the number of coding exons
	4. If there is a tie then look at the difference between the number of exons
	5. If there is a tie then look at the difference between the total coding length
	6. If still tied then pick all of them and report the situation
	"""
	
	oneToOneTrs=[]
	newTsimMatrixAll=tsimMatrixAll.copy()
	newTsimMatrixCoding=tsimMatrixCoding.copy()

	noTr1s=len(tdic1)
	noTr2s=len(tdic2)
	newMax=0.001
	conditionCounter=[0,0,0,0,0,0,0] #1=Unique winner, 2=tie in one score, not in the other, 
	#3= tie in both scores but coding exon length diff breaks the tie, 
	#4= tie in both scores and coding exon length diff but overall exon length breaks the tie
	#5= tie in all the above but coding length (bp) diff breaks the tie
	#6= tie in all the above, just give up and report all

	## I dropped the extra condition that ensured not too many mappings so that now ties are allowed to be reported
	#while sum(sum(newTsimMatrixCoding))>0 and newMax>0 and len(oneToOneTrs)<min(noTr1s,noTr2s):
	while sum(sum(newTsimMatrixAll))+sum(sum(newTsimMatrixCoding))>0 and newMax>0:


		# select the max among the two score matrices
		newMaxAll=newTsimMatrixAll.max()
		newMaxCoding=newTsimMatrixCoding.max()
		newMax=max(newMaxCoding,newMaxAll)
		

		# This is a threshold choice but it makes sense to break when transcript similarity scores are too low 
		if newMaxAll<0.2: 
			break

		codingIsMax=True 
		# choose where the max of the two came from. For precision purposes use a small difference instead of equality check.
		if newMaxAll>newMaxCoding:
			xs,ys=np.where(abs(newTsimMatrixAll-newMax)<=0.00000000001)
			codingIsMax=False
		else:
			xs,ys=np.where(abs(newTsimMatrixCoding-newMax)<=0.00000000001)
			codingIsMax=True
		#

		# Condition 0: something wrong, shouldn't have entered the while loop
		if len(xs)==0 or len(ys)==0: 
			conditionCounter[0]+=1

		# Condition 1: one unique maximum scoring transcript pair
		if len(xs)==1: 
			#print "Condition1"
			x=xs[0]; y=ys[0] 
			oneToOneTrs.append([x,y,1,newTsimMatrixCoding[x,y],newTsimMatrixAll[x,y],-1,-1])
			newTsimMatrixCoding[x,:]=0; newTsimMatrixCoding[:,y]=0
			newTsimMatrixAll[x,:]=0; newTsimMatrixAll[:,y]=0
			conditionCounter[1]+=1
# This automatically makes anything below unreachable and just reports the first out of ties #######################
#		elif len(xs)>1: 
#			x=xs[0]; y=ys[0] 
#			oneToOneTrs.append([x,y,1,newTsimMatrixCoding[x,y],newTsimMatrixAll[x,y],-1,-1])
#			newTsimMatrixCoding[x,y]=0; newTsimMatrixCoding[x,y]=0
#			newTsimMatrixAll[x,y]=0; newTsimMatrixAll[x,y]=0
#			conditionCounter[1]+=1
####################################################################################################################

		# Condition 1 fails -> multiple equal maximum similarity (either coding or all)
		else: 
			#print "Condition1+"
			allMax=0
			axs,ays=[],[]
			# loop and find the max among the xs,ys pairs
			for i in range(len(xs)):
				x=xs[i]; y=ys[i]
				if codingIsMax==True:
					if newTsimMatrixAll[x,y]>allMax:
						allMax=newTsimMatrixAll[x,y]
				else:
					if newTsimMatrixCoding[x,y]>allMax:
						allMax=newTsimMatrixCoding[x,y]
				#
			#

			# loop once more and find the pairs among xs,ys that still tie
			for i in range(len(xs)):
				x=xs[i]; y=ys[i]
				if codingIsMax==True:
					if abs(newTsimMatrixAll[x,y]-allMax)<=0.00000000001:
						axs.append(x); ays.append(y)
				else:
					if abs(newTsimMatrixCoding[x,y]-allMax)<=0.00000000001:
						axs.append(x); ays.append(y)
			#

			# Condition 2: one unique maximum in one matrix after tie in the other
			if len(axs)==1: 
				#print ("Condition2")
				x=axs[0]; y=ays[0]
				oneToOneTrs.append([x,y,2,newTsimMatrixCoding[x,y],newTsimMatrixAll[x,y],-1,-1])
				newTsimMatrixCoding[x,:]=0; newTsimMatrixCoding[:,y]=0
				newTsimMatrixAll[x,:]=0; newTsimMatrixAll[:,y]=0
				conditionCounter[2]+=1

			# Condition 2 fails -> multiple equal coding and overall sim. check the number of coding exons for tie break
			else: 
				#print ("Condition2+")
				minCodingExonDiffList=[]
				for i in range(len(axs)):
					x=axs[i]; y=ays[i]
					t1=trevDic1[x]; t2=trevDic2[y]; 
					tempCodingDiff=abs(len(transcriptDic1[t1].codingExons)-len(transcriptDic2[t2].codingExons))
					minCodingExonDiffList.append(tempCodingDiff)
				#
				minCodingExonDiffList=np.array(minCodingExonDiffList)
				minCodingExonDiff=np.min(minCodingExonDiffList)
				ncis=np.where(minCodingExonDiffList==minCodingExonDiff)[0] # one-sided index

				# Condition 3: one unique minimum difference in the number of coding exons
				if len(ncis)==1: 
					#print ("Condition3")
					x=axs[ncis[0]]; y=ays[ncis[0]] 
					oneToOneTrs.append([x,y,3,newTsimMatrixCoding[x,y],newTsimMatrixAll[x,y],minCodingExonDiff,-1])
					newTsimMatrixCoding[x,:]=0; newTsimMatrixCoding[:,y]=0
					newTsimMatrixAll[x,:]=0; newTsimMatrixAll[:,y]=0
					conditionCounter[3]+=1
				# Condition 3 fails -> multiple equal coding/ and overall sim and coding exon count difference
				else: 
					#print ("Condition3+")
					minAllExonDiffList=[]
					for i in ncis:
						x=axs[i]; y=ays[i]
						t1=trevDic1[x]; t2=trevDic2[y]
						tempAllDiff=abs(len(transcriptDic1[t1].exons)-len(transcriptDic2[t2].exons))
						minAllExonDiffList.append(tempAllDiff)
					#
					minAllExonDiffList=np.array(minAllExonDiffList)
					minAllExonDiff=np.min(minAllExonDiffList)
					nais=np.where(minAllExonDiffList==minAllExonDiff)[0] # one-sided index

					# Condition 4: one unique minimum difference in the number of all exons
					if len(nais)==1: 
						#print ("Condition4")
						x=axs[ncis[nais[0]]]; y=ays[ncis[nais[0]]]
						oneToOneTrs.append([x,y,4,newTsimMatrixCoding[x,y],newTsimMatrixAll[x,y],minCodingExonDiff,minAllExonDiff])
						newTsimMatrixCoding[x,:]=0; newTsimMatrixCoding[:,y]=0
						newTsimMatrixAll[x,:]=0; newTsimMatrixAll[:,y]=0
						conditionCounter[4]+=1
					# Condition 4 fails -> multiple equal coding/ and overall sim andcoding/ and overall exon count difference
					else: 
						#print ("Condition4+")
						# check the coding lengths of transcripts
						minCodingLenDiffList=[]
						for i in nais:
							x=axs[ncis[i]]; y=ays[ncis[i]]
							t1=trevDic1[x]; t2=trevDic2[y]; 
							tlen1=0; tlen2=0
							for stC,enC in transcriptDic1[t1].codingExons:
								tlen1+=abs(int(stC)-int(enC))
							for stC,enC in transcriptDic2[t2].codingExons:
								tlen2+=abs(int(stC)-int(enC))
							tempLenDiff=abs(tlen1-tlen2)
							minCodingLenDiffList.append(tempLenDiff)
						#
						minCodingLenDiffList=np.array(minCodingLenDiffList)
						minCodingLenDiff=np.min(minCodingLenDiffList)
						lcis=np.where(minCodingLenDiffList==minCodingLenDiff)[0] # one-sided index

						# Condition 5: one unique minimum difference in the coding length differences
						if len(lcis)==1: 
							#print ("Condition5")
							x=axs[ncis[nais[lcis[0]]]]; y=ays[ncis[nais[lcis[0]]]]
							oneToOneTrs.append([x,y,5,newTsimMatrixCoding[x,y],newTsimMatrixAll[x,y],minCodingExonDiff,minAllExonDiff])
							newTsimMatrixCoding[x,:]=0; newTsimMatrixCoding[:,y]=0
							newTsimMatrixAll[x,:]=0; newTsimMatrixAll[:,y]=0
							conditionCounter[5]+=1
						# Condition 5 fails -> all the above could not break the tie. Just report all of them.
						else: 
							#print (["Condition5+ = 6", lcis, len(lcis),len(oneToOneTrs)])
							for i in lcis:
								x=axs[ncis[nais[i]]]; y=ays[ncis[nais[i]]]
								oneToOneTrs.append([x,y,6,newTsimMatrixCoding[x,y],\
									newTsimMatrixAll[x,y],minCodingExonDiff,minAllExonDiff])
							for i in lcis:
								x=axs[ncis[nais[i]]]; y=ays[ncis[nais[i]]]
								newTsimMatrixCoding[x,:]=0; newTsimMatrixCoding[:,y]=0
								newTsimMatrixAll[x,:]=0; newTsimMatrixAll[:,y]=0
							#
							conditionCounter[6]+=1
					#
			#
		#
	#
	print ("\nCondition counter from the greedy transcript mapping stage:")
	print ("\t%d pairs with Condition1: Unique winner pair" % conditionCounter[1])
	print ("\t%d pairs with Condition2: Tie in one score, not in the other" % conditionCounter[2])
	print ("\t%d pairs with Condition3: Tie in both scores but coding exon length diff breaks the tie" % conditionCounter[3])
	print ("\t%d pairs with Condition4: Tie in both scores and coding exon length diff but overall exon length breaks the tie" % conditionCounter[4])
	print ("\t%d pairs with Condition5: Tie in all the above but coding length (bp) diff breaks the tie" % conditionCounter[5])
	print ("\t%d pairs with Condition6: Tie in all the above, just give up and report all" % conditionCounter[6])

	return oneToOneTrs

def greedy_exonPairs_per_transcriptPair(t1,t2):
	"""
	For a given transcript pair find a matching of exons greedily by picking the
	highest similarity exon pairs first and so on to compute one-to-one exon mappings 
	and an overall pairwise transcript similarity score.
	"""

	# list the exons and coding exons 
	elist1=transcriptDic1[t1].exons
	elist2=transcriptDic2[t2].exons

	#ecodinglist1=transcriptDic1[t1].codingExons # FERHAT fixed 4/5/2018
	ecodinglist1=[]
	for ex1 in [x[2] for x in elist1]:
		if exonDic1[ex1].exon_type!="nonCoding": ecodinglist1.append(ex1)

	#ecodinglist2=transcriptDic2[t2].codingExons # FERHAT fixed 4/5/2018
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

	# reverse dictionaries for lookup from index to exon id
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


	# greedy selection of best exon pairs using the two exon similarty matrices (coding and all)
	newEsimMatrixAll=esimMatrixAll.copy()
	newEsimMatrixCoding=esimMatrixCoding.copy()
	newMax=0.001
	conditionCounter=[0,0,0,0] #1=Unique winner, 2=breakable ties, 3=unbreakable ties
	while sum(sum(newEsimMatrixAll))+sum(sum(newEsimMatrixCoding))>0 and newMax>0 and len(oneToOneExs)<min(len(elist1),len(elist2)):
		# select the max among the two score matrices
		newMaxAll=newEsimMatrixAll.max()
		newMaxCoding=newEsimMatrixCoding.max()
		newMax=max(newMaxAll,newMaxCoding)

		codingIsMax=True 
		# choose where the max of the two came from. For precision purposes use a small difference instead of equality check.
		if newMaxAll>newMaxCoding:
			xs,ys=np.where(abs(newEsimMatrixAll-newMax)<=0.00000000001)
			codingIsMax=False
		else:
			xs,ys=np.where(abs(newEsimMatrixCoding-newMax)<=0.00000000001)
			codingIsMax=True

		# Condition 0: something wrong, shouldn't have entered the while loop
		if len(xs)==0 or len(ys)==0: 
			conditionCounter[0]+=1

		# Condition 1: one unique maximum for the maximum similarity (either coding or all)
		if len(xs)==1: 
			x=xs[0]; y=ys[0] 
			oneToOneExs.append([x,y,1,newEsimMatrixCoding[x,y],newEsimMatrixAll[x,y]])
			newEsimMatrixCoding[x,:]=0; newEsimMatrixCoding[:,y]=0
			newEsimMatrixAll[x,:]=0; newEsimMatrixAll[:,y]=0
			conditionCounter[1]+=1
		# Condition 1 fails -> multiple equal maximum similarity (either coding or all)
		else: 
			allMax=0
			axs,ays=[],[]
			# loop and find the max among the xs,ys pairs
			for i in range(len(xs)):
				x=xs[i]; y=ys[i]
				if codingIsMax==True:
					if newEsimMatrixAll[x,y]>allMax:
						allMax=newEsimMatrixAll[x,y]
				else:
					if newEsimMatrixCoding[x,y]>allMax:
						allMax=newEsimMatrixCoding[x,y]
				#
			#
			# loop once more and find the pairs among xs,ys that still tie
			for i in range(len(xs)):
				x=xs[i]; y=ys[i]
				if codingIsMax==True:
					if abs(newEsimMatrixAll[x,y]-allMax)<=0.00000000001:
						axs.append(x); ays.append(y)
				else:
					if abs(newEsimMatrixCoding[x,y]-allMax)<=0.00000000001:
						axs.append(x); ays.append(y)
			#

			# Condition 2: one unique maximum in 'the other' matrix after tie in max (either coding or all)
			if len(axs)==1: 
				x=axs[0]; y=ays[0] 
				oneToOneExs.append([x,y,2,newEsimMatrixCoding[x,y],newEsimMatrixAll[x,y]])
				newEsimMatrixCoding[x,:]=0; newEsimMatrixCoding[:,y]=0
				newEsimMatrixAll[x,:]=0; newEsimMatrixAll[:,y]=0
				conditionCounter[2]+=1
			# Condition 2 fails ->  multiple equal coding and overall sim - THEN JUST pick the first one and report it
			else: 
				# Condition 3: JUST pick the first one and report it
				x=axs[0]; y=ays[0] 
				oneToOneExs.append([x,y,3,newEsimMatrixCoding[x,y],newEsimMatrixAll[x,y]])
				newEsimMatrixCoding[x,:]=0; newEsimMatrixCoding[:,y]=0
				newEsimMatrixAll[x,:]=0; newEsimMatrixAll[:,y]=0
				conditionCounter[3]+=1
		#
	#

	exscoreAll=0
	exscoreCoding=0
	expathAll=[]
	expathCoding=[]
	for l in oneToOneExs:
		ei,ej,condition,codingScore,allScore=l
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

	## taking max of exon lengths
	#exscoreAll=(exscoreAll)/max(len(elist1),len(elist2))
	#exscoreCoding=(exscoreCoding)/max(1,max(len(ecodinglist1),len(ecodinglist2)))

	## taking sum of exon lengths
	exscoreAll=(2*exscoreAll)/max(1,len(elist1)+len(elist2))
	exscoreCoding=(2*exscoreCoding)/max(1,len(ecodinglist1)+len(ecodinglist2))

	return (exscoreAll,expathAll,exscoreCoding,expathCoding)


def output_transcript_level_matches(genePairId,exonClasses1,exonClasses2):
	"""
	Find exon level pairings/mappings using previously calculated exon similarities
	and a method of choice (greedy/DP/allToAll) to compute transcript similarities.
	Then, extract transcript level mappings using a greedy algorithm that also
	does tie breaking in a rule-based system. 

	"""	

	print ("*****************************************************************")
	outfilename=outdir+"/exonLevelMappings-"+str(MAPPED_EXON_THRES)+".txt"
	print ("Writing exon-level mappings into file:\n %s" % outfilename)
	outfileExonMappings=open(outfilename,'w')
	outfileExonMappings.write("chrName1\tstartCoord1\tendCoord1\tstrand1\texonID1\texonType1\tchrName2\tstartCoord2\tendCoord2\tstrand2\texonID2\texonType2\toverlapScoreFromFullLength\toverlapScoreFromPartialCodingPart\n")

	outfilename=outdir+"/transcriptLevelSimilarities-"+str(MAPPED_EXON_THRES)+".txt"
	print ("Writing trascript-level similarity scores into file:\n %s" % outfilename)
	outfileTranscriptSims=open(outfilename,'w')
	outfileTranscriptSims.write("chrName1\tstartCoord1\tendCoord1\tstrand1\ttranscriptID1\ttranscriptType1\tchrName2\tstartCoord2\tendCoord2\tstrand2\ttranscriptID2\ttranscriptType2\tnoOfAllExons1\tnoOfAllExons2\tmappedAllExonPairs\tnoOfCodingExons1\tnoOfCodingExons2\tmappedCodingExonPairs\toverallSimScore\tcodingSimScore\n")


	# list of all transcripts per each gene
	tdic1=geneDic1[g1].transcripts
	tdic2=geneDic2[g2].transcripts
	
	# pairwise transcript similarity matrices
	noRows=len(tdic1)
	noCols=len(tdic2)
	tsimMatrixAll=np.zeros((noRows,noCols))
	tsimMatrixCoding=np.zeros((noRows,noCols))

	allExonPairs={}
	ti=0
	for t1 in tdic1:
		elist1=transcriptDic1[t1].exons
		#ecodinglist1=transcriptDic1[t1].codingExons # FERHAT fixed 4/5/2018
		ecodinglist1=[]
		for ex1 in [x[2] for x in elist1]:
			if exonDic1[ex1].exon_type!="nonCoding": ecodinglist1.append(ex1)


		t1Entry=transcriptDic1[t1]
		ch1,s1="chr"+t1Entry.basicInfoDic["chromosome"],t1Entry.basicInfoDic["start_coord"]
		e1,strand1=t1Entry.basicInfoDic["end_coord"],t1Entry.basicInfoDic["strand"]
		t1type=t1Entry.basicInfoDic["transcript_biotype"]

		tj=0
		for t2 in tdic2:
			elist2=transcriptDic2[t2].exons
			#ecodinglist2=transcriptDic2[t2].codingExons # FERHAT fixed 4/5/2018
			ecodinglist2=[]
			for ex2 in [x[2] for x in elist2]:
				if exonDic2[ex2].exon_type!="nonCoding": ecodinglist2.append(ex2)

			t2Entry=transcriptDic2[t2]
			ch2,s2="chr"+t2Entry.basicInfoDic["chromosome"],t2Entry.basicInfoDic["start_coord"]
			e2,strand2=t2Entry.basicInfoDic["end_coord"],t2Entry.basicInfoDic["strand"]
			t2type=t2Entry.basicInfoDic["transcript_biotype"]

			### extract exon mappings/pairings using one of the below methods to compute transcript level similarities
			allScore,pathAll,codingScore,pathCoding=greedy_exonPairs_per_transcriptPair(t1,t2)
			#allScore,pathAll,codingScore,pathCoding=dynamic_programming_per_transcriptPair(t1,t2)
			#allScore,pathAll,codingScore,pathCoding=all_exonPairs_per_transcriptPair(t1,t2)

			outfileTranscriptSims.write("%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t\t%.3f\t%.3f\n" % \
				(ch1,s1,e1,strand1,t1,t1type,ch2,s2,e2,strand2,t2,t2type, \
				len(elist1),len(elist2),len(pathAll), len(ecodinglist1),len(ecodinglist2), len(pathCoding), \
				allScore,codingScore))

			tsimMatrixAll[ti,tj]=allScore
			tsimMatrixCoding[ti,tj]=codingScore

			for l in pathAll:
				ei,ej,allOverlapScore,codingOverlapScore=l[0],l[1],l[2],l[3]
				ex1=elist1[ei][2]
				ex2=elist2[ej][2]
				if ex1 not in allExonPairs:
					allExonPairs[ex1]={}
				if ex2 not in allExonPairs[ex1]:
					allExonPairs[ex1][ex2]=[0,0]
				allExonPairs[ex1][ex2]=[allOverlapScore,codingOverlapScore]
			#
			tj+=1

		#
		ti+=1
	#
	for ex1 in allExonPairs:
		origExonEntry1=exonDic1[ex1]
		ch1,s1="chr"+origExonEntry1.basicInfoDic["chromosome"],origExonEntry1.basicInfoDic["start_coord"]
		e1,strand1=origExonEntry1.basicInfoDic["end_coord"],origExonEntry1.basicInfoDic["strand"]
		e1type=origExonEntry1.exon_type
		for ex2 in allExonPairs[ex1]:
			origExonEntry2=exonDic2[ex2]
			ch2,s2="chr"+origExonEntry2.basicInfoDic["chromosome"],origExonEntry2.basicInfoDic["start_coord"]
			e2,strand2=origExonEntry2.basicInfoDic["end_coord"],origExonEntry2.basicInfoDic["strand"]
			e2type=origExonEntry2.exon_type
			allOverlapScore,codingOverlapScore=allExonPairs[ex1][ex2]
			# Only output exon mappings that pass the threshold even though other pairs are also used for
			# computing the transcript level similarities
			if allOverlapScore>=MAPPED_EXON_THRES or codingOverlapScore>=MAPPED_EXON_THRES:
				outfileExonMappings.write("%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%.3f\t%.3f\n" \
					% (ch1,s1,e1,strand1,ex1,e1type,ch2,s2,e2,strand2,ex2,e2type,allOverlapScore,codingOverlapScore))
		#
	#
	outfileExonMappings.close()
	outfileTranscriptSims.close()

	outfilename=outdir+"/transcriptLevelMappings-"+str(MAPPED_EXON_THRES)+".txt"
	print ("Writing transcript-level mappings into file:\n %s" % outfilename)
	outfileTranscriptMappings=open(outfilename,'w')
	outfileTranscriptMappings.write("chrName1\tstartCoord1\tendCoord1\tstrand1\ttranscriptID1\ttranscriptType1\tchrName2\tstartCoord2\tendCoord2\tstrand2\ttranscriptID2\ttranscriptType2\tnoOfAllExons1\tnoOfAllExons2\tmappedAllExonPairs\tnoOfCodingExons1\tnoOfCodingExons2\tmappedCodingExonPairs\toverallSimScore\tcodingSimScore\n")


	trevDic1={}
	c=0
	for t in tdic1:
		trevDic1[c]=t
		c+=1 
	c=0
	trevDic2={}
	for t in tdic2:
		trevDic2[c]=t
		c+=1 

	# find and write the one-to-one mapped transcript pairs
	oneToOneTrs=greedily_map_transcripts_with_tie_breaks(tdic1,tdic2,tsimMatrixAll,tsimMatrixCoding,trevDic1,trevDic2)
	for l in oneToOneTrs:
		ti,tj,condition,codingScore,allScore,codingExonDiff,allExonDiff=l
		t1=trevDic1[ti]
		t2=trevDic2[tj]

		elist1=transcriptDic1[t1].exons
		#ecodinglist1=transcriptDic1[t1].codingExons # FERHAT fixed 4/5/2018
		ecodinglist1=[]
		for ex1 in [x[2] for x in elist1]:
			if exonDic1[ex1].exon_type!="nonCoding": ecodinglist1.append(ex1)

		t1Entry=transcriptDic1[t1]
		ch1,s1="chr"+t1Entry.basicInfoDic["chromosome"],t1Entry.basicInfoDic["start_coord"]
		e1,strand1=t1Entry.basicInfoDic["end_coord"],t1Entry.basicInfoDic["strand"]
		t1type=t1Entry.basicInfoDic["transcript_biotype"]

		elist2=transcriptDic2[t2].exons
		#ecodinglist2=transcriptDic2[t2].codingExons # FERHAT fixed 4/5/2018
		ecodinglist2=[]
		for ex2 in [x[2] for x in elist2]:
				if exonDic2[ex2].exon_type!="nonCoding": ecodinglist2.append(ex2)

		t2Entry=transcriptDic2[t2]
		ch2,s2="chr"+t2Entry.basicInfoDic["chromosome"],t2Entry.basicInfoDic["start_coord"]
		e2,strand2=t2Entry.basicInfoDic["end_coord"],t2Entry.basicInfoDic["strand"]
		t2type=t2Entry.basicInfoDic["transcript_biotype"]

		outfileTranscriptMappings.write("%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t\t%.3f\t%.3f\n" % \
			(ch1,s1,e1,strand1,t1,t1type,ch2,s2,e2,strand2,t2,t2type, \
			len(elist1),len(elist2),len(pathAll), len(ecodinglist1),len(ecodinglist2), len(pathCoding), \
			allScore,codingScore))
	#
	outfileTranscriptMappings.close()

	#Write bed files per pair to be able to visualize in genome browser
	write_out_bedFiles_for_UCSCbrowser(genePairId,exonClasses1,exonClasses2,oneToOneTrs,trevDic1,trevDic2)
	return

def classify_exons(genePairId,mappedExonThreshold):
	"""
	Given a threshold on the similarity between two exons to be
	deemed "mapped" (mappedExonThreshold), this function determines
	the exon classes using all the info from parsed liftover files.
	"""
#	print ("*****************************************************************")
	# open the output files to which the results of the mapping will be written
	outfilename=outdir+"/exonLevelSimilarities-"+str(MAPPED_EXON_THRES)+".txt"
	outfileExonSims=open(outfilename,'w')
	print ("Writing exon-level similarity scores into file:\n %s" % outfilename)

	outfileExonSims.write("chrName1\tstartCoord1\tendCoord1\tstrand1\texonID1\texonType1\tchrName2\tstartCoord2\tendCoord2\tstrand2\texonID2\texonType2\toverlapScoreFromFullLength\toverlapScoreFromPartialCodingPart\n")

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
				outfileExonSims.write("%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%.3f\t%.3f\n" \
					% (ch1,s1,e1,strand1,ex1,e1type,ch2,s2,e2,strand2,ex2,e2type,allOverlapScore,codingOverlapScore)),
			#
			#outfileExonSims.write("%s\t%s\t%s\t%s\t%.3f\t%.3f\t%s\t%s\n" % \
			#	(ex1,ex2,origExonEntry1.get_summary_string(),origExonEntry2.get_summary_string(),allOverlapScore,\
			#	codingOverlapScore,minMatchStrAll,minMatchStrCoding))
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
			#outfileExonSims.write("PotentialSpliceLoss\t%r\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%d\n" % \
			#	(junctionProblem,minMatchStr,origExonEntry.get_summary_string(),ch,int(s),int(e),strand,sqAcc,sqDon,liftedOverLenForEx))
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
			#outfileExonSims.write("UnmappedExon\t%s\t%s\t%s\n" % (minMatchStr,origExonEntry.get_summary_string(),whyUnmapped))
		#
	#
	outfileExonSims.close()


	outfilename=outdir+"/exonClasses-"+str(MAPPED_EXON_THRES)+".txt"
	outfileExonClasses=open(outfilename,'w')
	print ("\nWriting exon classes into file:\n %s" % outfilename)

	outfileExonClasses.write("exonID\texonClass\tmaxOverlapScore\tmaxCodingOverlapScore\twhyUnmapped\tjunctionProblem\tchrName\tstartCoord\tendCoord\tstrand\texonID\texonType\tcodingStart\tcodingEnd\ttranscriptIDs\texonNumbers\tgeneID\texonLength\tacceptor2bp\tdonor2bp\tavgCodingConsScore\tavgConsScore\tfirstMidLastCounts\telementType\n")

	liftOverExpDic={'#Partially': "PartiallyDeletedInNew", "#Split": "SplitInNew", "#Deleted" : "DeletedInNew", \
		"#Duplicated": "DuplicatedInNew", "#Boundary" : "BoundaryProblem"}	

	#print("***************** Exon classes for the first organism **************************\n"),
	counts=[0,0,0,0]
	for ex in exonDic1:
		origExonEntry=exonDic1[ex]
		etype=origExonEntry.exon_type
		allOverlapScore=0.0; codingOverlapScore=0.0; whyUnmapped='NA'; junctionProblem='NA'
		cl="OTHER"
		if ex in exonClasses1:
			cl,field1,field2=exonClasses1[ex][0],exonClasses1[ex][1],exonClasses1[ex][2]
			if cl=='Mapped': 
				counts[0]+=1
				allOverlapScore=float(field1); codingOverlapScore=float(field2)
				if etype=='fullCoding': 
					codingOverlapScore=max(allOverlapScore,codingOverlapScore)
					allOverlapScore=codingOverlapScore
			elif cl=='Unmapped':
				counts[1]+=1
				whyUnmapped=liftOverExpDic[field1]
			elif cl=='Nonintersecing':
				counts[2]+=1
				junctionProblem=field1
		#
		if cl=='OTHER': 
			counts[3]+=1
		outfileExonClasses.write("%s\t%s\t%.3f\t%.3f\t%s\t%s\t%s\n" % \
			(ex,cl,allOverlapScore,codingOverlapScore,whyUnmapped,junctionProblem,origExonEntry.get_summary_string()))
	#
	print ("\tFor org1: Mapped exons= %d, Unmapped exons= %d, Nonintersecting exons= %d, OTHER= %d" % (counts[0],counts[1],counts[2],counts[3]))

	#print("***************** Exon classes for the second organism **************************\n"),
	counts=[0,0,0,0]
	for ex in exonDic2:
		origExonEntry=exonDic2[ex]
		etype=origExonEntry.exon_type
		allOverlapScore=0.0; codingOverlapScore=0.0; whyUnmapped='NA'; junctionProblem='NA'

		cl="OTHER"
		if ex in exonClasses2:
			cl,field1,field2=exonClasses2[ex][0],exonClasses2[ex][1],exonClasses2[ex][2]
			if cl=='Mapped': 
				counts[0]+=1
				allOverlapScore=float(field1); codingOverlapScore=float(field2)
				if etype=='fullCoding': 
					codingOverlapScore=max(allOverlapScore,codingOverlapScore)
					allOverlapScore=codingOverlapScore
			elif cl=='Unmapped':
				counts[1]+=1
				whyUnmapped=liftOverExpDic[field1]
			elif cl=='Nonintersecing':
				counts[2]+=1
				junctionProblem=field1
		#
		if cl=='OTHER': 
			counts[3]+=1
		outfileExonClasses.write("%s\t%s\t%.3f\t%.3f\t%s\t%s\t%s\n" % \
			(ex,cl,allOverlapScore,codingOverlapScore,whyUnmapped,junctionProblem,origExonEntry.get_summary_string()))
	#	
	print ("\tFor org2: Mapped exons= %d, Unmapped exons= %d, Nonintersecting exons= %d, OTHER= %d" % (counts[0],counts[1],counts[2],counts[3]))
	print ("*****************************************************************\n")
	#
	outfileExonClasses.close()
	
	return exonClasses1,exonClasses2

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
				continue
			elif max(0, min(int(e),int(gE))-max(int(s),int(gS)))<1:
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

			## below check not needed anymore due to above ch==gCh check
			#if ch not in genome:
			#	continue

			sq=genome[ch].seq[int(s)-3:int(e)+2].tostring().lower().upper().decode("utf-8") 
			print (sq)
			if strand=="-":
				sqAcc=sq[-2:].lower().translate(complement)[::-1].upper() # reverse complement
				sqDon=sq[:2].lower().translate(complement)[::-1].upper() # reverse complement
			else:
				sqAcc=sq[:2]
				sqDon=sq[-2:]
			#

			# Below if else structure is used for determining junction problems 
			# depending on whether an exon is first, mid, last or a combination 
			f,m,l=fml1
			junctionProblem=False
			if m>0 or (f>0 and l>0):
				if sqAcc!=acceptor1 or sqDon!=donor1:
					junctionProblem=True
			elif f>0:
				if sqDon!=donor1:
					junctionProblem=True
			elif l>0:
				if sqAcc!=acceptor1:
					junctionProblem=True
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
			#print ("PotentialSpliceLoss\t%r\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%d\n" % \
			#	(junctionProblem,minMatchStr,origExonEntry.get_summary_string(),ch,int(s),int(e),strand,sqAcc,sqDon,liftedOverLenForEx)),
		#
	#
	infile.close()
	return

def parse_liftOver_unmapped_file_perExon(infilename):
	"""
	Parse the .unmapped files.
	An example line:  
		chrX	100593624	100594035	+	ENSE00001952391	#Split	flank0-minMatch1.0-multiples
	# no problem with symmetry, just add all unmapped exons from org1 and org2 into unmappedExons dictionary
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
		#print ("UnmappedExon\t%s\t%s\t%d\t%d\t%s\t%s\t%s\n" % (ex,ch,int(s),int(e),strand,whyUnmapped,minMatchStr)),
	#
	infile.close()
	return

def compute_overlap_score_for_exons(exonLength1, exonLength2, liftedOverLength1, overlap):
	r=1 # if exonLength1==liftedOverLength1, meaning liftOver did not change the overall length
	#print [overlap, liftedOverLength1, exonLength2, exonLength1]
	overlapScore=(2.0*overlap)/(liftedOverLength1+exonLength2)
	if exonLength1<liftedOverLength1:
		r=(1.0*exonLength1)/liftedOverLength1
	elif exonLength1>liftedOverLength1:
		r=(1.0*liftedOverLength1)/exonLength1
	return min(1.0, r*overlapScore)

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
		else: # entry is reporting overlap between two full length exons (either fullCoding or nonCoding or partCoding)
			allOverlapScore=compute_overlap_score_for_exons(origExonLen1,origExonLen2,liftedOverLenForEx,overlap)
		#
		# make sure to set coding overlap to the same number as overall overlap for full coding pairs
		# if origExonEntry1.exon_type=="fullCoding" and origExonEntry2.exon_type=="fullCoding":
		# codingOverlapScore=allOverlapScore

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

def parse_liftOver_all_files_perGenePair(genePairId):
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
	print ("*****************************************************************")
	for exonDic in [exonDic1,exonDic2]:
		isOneToTwo = True if c==1 else False #ifelse one liner
		counts={"000": 0, "001": 0, "010": 0, "011": 0, "100": 0, "101": 0, "110": 0, "111": 0}
		#
		for exonId in exonDic:
			infilenameMapped=perExonLiftoverDir+"/org"+str(c)+"/"+exonId+"_mapped.txt"
			infilenameUnmapped=perExonLiftoverDir+"/org"+str(c)+"/"+exonId+"_unmapped.txt"
			infilenameNonintersecting=perExonLiftoverDir+"/org"+str(c)+"/"+exonId+"_nonintersecting.txt"

			isM,isU,isN=os.path.isfile(infilenameMapped),os.path.isfile(infilenameUnmapped),os.path.isfile(infilenameNonintersecting)
			#       counts["001"]=counts["001"]+1 # only nonintersecting - ~6k
			#       counts["010"]=counts["010"]+1 # only unmapped - ~69k
			#       counts["011"]=counts["011"]+1 # unmapped and nonintersecting - ~13k
			#       counts["100"]=counts["100"]+1  # only mapped - ~258k
			#       counts["101"]=counts["101"]+1 # mapped and nonintersecting - ~60k
			#       counts["110"]=counts["110"]+1 # mapped and unmapped - ~83k
			#       counts["111"]=counts["111"]+1 # all three files - ~11k

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
		if c==1: print ("Exon file type summaries for the first gene from: "+str(genePairId))
		else: print ("Exon file type summaries for the second gene from: "+str(genePairId))
		#k=["000", "001", "010", "011", "100", "101", "110", "111"]
		#desc=["No file exists", "Only nonintersecting", "Only unmapped", "nonintersecting and unmapped",\
		#"Only Mapped", "mapped and nonintersecting", "mapped and unmapped", "All three files"]
		k=[ "000", "100", "001", "010", "110", "101", "011", "111"]
		desc=["No file exists", "Only Mapped", "Only nonintersecting", "Only unmapped", \
		"Mapped and unmapped", "Mapped and nonintersecting", "Nonintersecting and unmapped", "All three files"]

		for x in range(len(k)):
			print ("\t%d exons with: %s" % (counts[k[x]],desc[x]))
		print
		c+=1
	#
	return # return from parse_liftOver_all_files_perGenePair

def sort_exonDics_byCoordinates(origExonCoords):
	"""
	Order the exons wrt their start coordinates and delete the ones
	that are duplicates in terms of coordinates. Also return a dict of
	duplicates.
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


def find_matchings_perGenePair(genePairId):
	"""
	This is the main function that takes a gene pair id and carries out
	all the steps for finding exon and transcript level mappings.
	"""
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
        ### exonDic1,exonDic2 have exons before duplicate removal in no specific order          ###
        ###########################################################################################

	print ("*****************************************************************")
	g1,g2=genePairId.split("-")
	print ("Gene pair ID: %s\n" % genePairId)
	print ("Information about each gene. Last two numbers are no of transcripts and exons")
	print ("%s\t%s\n" % (g1,geneDic1[g1].get_summary_string()))
	print ("%s\t%s\n" % (g2,geneDic2[g2].get_summary_string()))
	print ("Number of exons before and after duplicate removal according to coordinates")
	print ("Org1\t%d\t%d\n" % (len(exonDic1),len(origExonCoords1)))
	print ("Org2\t%d\t%d\n" % (len(exonDic2),len(origExonCoords2)))
	print ("*****************************************************************\n")

        # After this function call I will have all the information about mapped, unmapped, nonintersecting exons
        # in these three dictionaries: liftoverExonMappings, unmappedExons, nonintersectingExons
        # NOTE THAT these info will be at the exon level before duplicate removal!!!
	parse_liftOver_all_files_perGenePair(genePairId)

        # Now decide what class each exon will be assigned to among:
        #  Mapped, Unmapped, Nonintersecing or left unassigned (not in exonClasses)
	exonClasses1,exonClasses2=classify_exons(genePairId,MAPPED_EXON_THRES)

        #  Output transccript level matches in appropriate files
	output_transcript_level_matches(genePairId,exonClasses1,exonClasses2)

	return #exonLevelMatchScoresWithDups,transcriptLevelMatches

############## Function ends ##############

#   Testing functionalities
def main():

	###################################################################################
	## This is a parameter choice to determine what similarity level                 ##
	## (between 0 and 1) should be the threshold for mapping two exons to each other ##
	args = parse_args(sys.argv[1:])

	global MAPPED_EXON_THRES
	#MAPPED_EXON_THRES=float(argv[1])
	MAPPED_EXON_THRES=float(args.mapping)
	###################################################################################

	global indir,orgId1,orgId2,outdir
	# human interpretable names of input organisms that match those in config.conf file (e.g. human mouse)
	#orgId1=argv[2]
	#orgId2=argv[3]
	orgId1=args.org1
	orgId2=args.org2
	indir=ExTraMapperPath+"/preprocess/data/"+orgId1+"-"+orgId2

	# input directories that are assumed to exists
	global ensemblDir,genomedataDir
	global GTFsummaryDir,perGenePairPickleDir,liftOverFilesDir,perExonLiftoverDir
	ensemblDir=indir+"/ensemblDownloads"
	genomedataDir=indir+"/genomedataArchives"
	GTFsummaryDir=indir+"/GTFsummaries"
	perGenePairPickleDir=indir+"/perGenePairPickledInfo" # directory with one picked filed for gene pair information dictionaries
	liftOverFilesDir=indir+"/liftoverRelatedFiles"
	perExonLiftoverDir=indir+"/perExonLiftoverCoords" # directory for one file per exon liftover mappings
	#################

        # determine whether the user wants to run mapping for all possible pairs
        # or for just one pair that is given as input
	#if argv[4]=="all":
	#	genePairIds=sorted(os.listdir(perGenePairPickleDir))
	#else:
	#	genePairIds=[argv[4]]
	if args.ortholog == "all":
		genePairIds=sorted(os.listdir(perGenePairPickleDir))
	else:
		genePairIds=[args.ortholog]
        #

        ###################################################################################
        ### Some examples that can be inputted as individual pairs:
        ## In order to run one such example either do:
        # python ${EXTRAMAPPER_DIR}/scripts/ExTraMapper.py 0.9 $org1 $org2 ENSG00000000003-ENSMUSG00000067377
        ## OR simply remove the comment in front of it therefore overwriting what is given as argument
        ###################################################################################

        #genePairIds=["ENSG00000000003-ENSMUSG00000067377" # first gene pair in alphabetical order
        #genePairIds=["ENSG00000197410-ENSMUSG00000102692" # example of a chimeric mouse transcript that maps to two separate human transcripts
        #genePairIds=["ENSG00000000457-ENSMUSG00000026584" # example that I tested the duplicate removal part on
        #genePairIds=["ENSG00000124203-ENSMUSG00000050600" # one-to-one transcript matching with 5 exons on each side
        #genePairIds=["ENSG00000130368-ENSMUSG00000068037"] # this is 1 to 20 exon example that used to crash the code since adjacency matrices were all zeros
        #genePairIds=["ENSG00000121743-ENSMUSG00000048582"] # this is a 2 to 2 exon example that used to crash the code since simMatrix was all zeros
        #genePairIds=["ENSG00000197410-ENSMUSG00000102692"] # example of a chimeric mouse transcript that maps to two separate human transcripts
        #genePairIds=["ENSG00000146085-ENSMUSG00000023921"] # one-to-one transcript matching with 13 exons on each side
        #genePairIds=["ENSG00000185313-ENSMUSG00000034533"] # one-to-one transcript matching with 27 exons on one side 28 on the other
        #genePairIds=["ENSG00000004948-ENSMUSG00000023964"] # CALCR gene with deleted exons from human to mouse
        #genePairIds=["ENSG00000173769-ENSMUSG00000094985"] # TOPAZ1 gene with 20 exons and 1 transcript in each. 2nd exons say deleted but they blast.
        #genePairIds=["ENSG00000093167-ENSMUSG00000032497"] # an exon loss event in mouse due to splice junctions
        #genePairIds=["ENSG00000140416-ENSMUSG00000032366"] # TPM1 vs Tpm1
        #genePairIds=["ENSG00000078900-ENSMUSG00000029026"] # TP73
        #genePairIds=["ENSG00000049768-ENSMUSG00000039521"] # FOXP3
        #genePairIds=["ENSG00000157764-ENSMUSG00000002413"] # BRAF
        #genePairIds=["ENSG00000143315-ENSMUSG00000050229"] # 1 exon gene pair
        #genePairIds=["ENSG00000141510-ENSMUSG00000059552"] # TP53

        ## FIXME:delete below examples once done testing
        #genePairIds=["ENSG00000146085-ENSMUSG00000023921" # one-to-one transcript matching with 13 exons on each side
        ## THIS IS THE EXAMPLE that pointed out the importance of directions in the adjacency graphs ##
        ###############################################################################################
        ## THIS IS THE EXAMPLE that pointed out that some very small scores included heuristically in the matchings is not good
        #genePairIds="ENSG00000258052-ENSMUSG00000064181"] # one-to-one transcript matching with 10 exons on one side 5 on the other
        ###############################################################################################

       # This will run mappings in the all mode but only from a given start index to the end index
	#if argv[4]=="all" and len(argv)>5:
	#	startIndex=int(argv[5])
	#	endIndex=int(argv[6])
	#else: # either all or just one gene pair
	#	startIndex=0
	#	endIndex=len(genePairIds)
	startIndex=0
	endIndex=len(genePairIds)
        #

	c=0
	for genePairId in genePairIds[startIndex:endIndex]:
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

                # use only the ortho pairs with both sides protein_coding
		if origGeneEntry1.basicInfoDic["gene_biotype"]=="protein_coding" and origGeneEntry2.basicInfoDic["gene_biotype"]=="protein_coding":
			find_matchings_perGenePair(genePairId)
			c+=1
		else:
			sys.stderr.write("WARNING! Not running ExTraMapper for pair %s because one or both genes are not protein_coding\n" % (genePairId))
        #
	sys.stderr.write("\n........\nExTraMapper ran successfully for %d gene pairs between: %s and %s\n\n" % (c,orgId1,orgId2))

	return # from main


def parse_args(args):
	
	############################ Parse Arguments ############################
	parser = argparse.ArgumentParser(description="Check the help flag")
	parser.add_argument("-m", dest="mapping", help="ExTraMapper Exon threshold value [e.g. 1]", required=True)
	parser.add_argument("-o1", dest="org1", help="First organism name [e.g. human]", required=True)
	parser.add_argument("-o2", dest="org2", help="Second organism name [e.g. mouse]", required=True)
	parser.add_argument("-p", dest="ortholog", help="Orthologous gene pair [e.g. ENSG00000141510-ENSMUSG00000059552 OR all]", required=True)
	return parser.parse_args()

if __name__ == "__main__":
	main()
