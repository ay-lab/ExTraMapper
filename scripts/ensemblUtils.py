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
##############################################################
# Genomedata is off by one, seq[0] returns you the 1st bp. 
# Have to account for this by subtracting one from each seq coordinate
# before retrieving sequence from Genomedata.
from genomedata import Genome
##############################################################
import string
import math
import gzip
import cPickle as pickle
import numpy as np
complement = string.maketrans('atcgn', 'tagcn')

def parse_ensembl_geneAndProtein_pairings(infilename,proteinToGeneDic,proteinPairsDic):
	"""
	This function parses a given Ensembl file downloaded from below folders
	and parses out the gene and protein pairings. 
	For the first files, say org1 to org2, it only reports protein pairings.
	For the second file it reports back the protein and gene pairings using
	the protein pairings from the first proteinPairingsSoFar dictionary.
	Fields of these files are (as far as I understand):
		someScore chr	notsure	start end orthologytype
	
		0.14740	MT	192515	3307	4262	ortholog_one2one	73.39810	Euarchontoglires	77	
		ENSG00000198888	ENSP00000354687	ENSMUSP00000080991	77	0

	"""
	sys.stderr.write("Parsing gene and protein pairings from file "+infilename+"\n")

	isFirstFile=True
	if bool(proteinToGeneDic): # empty dic evaluates to False
		isFirstFile=False

	genePairsDic={}
	if infilename.endswith(".gz"):
		infile=gzip.open(infilename,'rb')
	else:
		infile=open(infilename,'r')
	#

	# it doesn't have a header
	lineCount=0
	for line in infile:
		words=line.rstrip().split()
		someScore,chr,notsure,st,en,orthologyType,someOtherScore,phylo,\
			gene1PercentIdentity,gene1,protein1,protein2,gene2PercentIdentity,notsure=words
		
		# skip the chromosomes that are not 1..23 or X or Y
		if chr=="MT" or len(chr)>2: 
			continue
		
		if "ortholog" not in orthologyType:
			continue

		proteinToGeneDic[protein1]=gene1
		if protein1 not in proteinPairsDic:
			proteinPairsDic[protein1]=[]
		proteinPairsDic[protein1].append(protein2)

		if not isFirstFile:
			if protein2 not in proteinToGeneDic: # second gene is not from a 1,22 or X,Y chr
				#print [chr,gene1,protein1,protein2]
				continue 
			gene2=proteinToGeneDic[protein2]
			# below checks ensure I only get one entry per g1-g2 pair
			if gene1 not in genePairsDic:
				genePairsDic[gene1]=[]
				genePairsDic[gene1].append([gene2,orthologyType,gene1PercentIdentity,gene2PercentIdentity])
			else: 
				if gene2 not in [a[0] for a in genePairsDic[gene1]]:
					genePairsDic[gene1].append([gene2,orthologyType,gene1PercentIdentity,gene2PercentIdentity])
			#
	#
	types={}
	for g1 in genePairsDic:
		type=genePairsDic[g1][0][1]
		num=len(genePairsDic[g1])
		if type=="ortholog_one2one" and num>1:
			sys.exit("Matching should be one2one but it isn't\t"+g1)
		if type not in types:
			types[type]=0
		types[type]+=1
	#
	infile.close()
	
	sys.stderr.write("Gene pair mappings summary: "+repr(types)+"\n\n")

	return proteinToGeneDic,genePairsDic,proteinPairsDic

def parse_ensembl_gene_pairings(infilename):
	"""
	This function parses a given Ensembl file (comma seperated) 
	that matches the genes of the first organism to the second.
	"""
	sys.stderr.write("Parsing gene pairings from file "+infilename+"\n")

	genePairsDic={}
	if infilename.endswith(".gz"):
		infile=gzip.open(infilename,'rb')
	else:
		infile=open(infilename,'r')
	#
	columnIndices={"GeneID1" : -1, "GeneID2" : -1, "Homology Type" : -1, "Orthology confidence": -1, "Percent Identity" : -1, "Chromosome Name" : []}
	# parse the information header first
	columnNames=infile.readline().strip().split(",")
	i=0 #0-based column indices
	for c in columnNames:
		if c=="Ensembl Gene ID":
			columnIndices["GeneID1"]=i
		elif c.endswith("Ensembl Gene ID"):
			columnIndices["GeneID2"]=i
		elif c=="Homology Type":
			columnIndices["Homology Type"]=i
		elif "Orthology confidence" in c:
			columnIndices["Orthology confidence"]=i
		elif c.endswith("Identity with respect to query gene"):
			columnIndices["Percent Identity"]=i
		elif c.endswith("Chromosome Name"):
			columnIndices["Chromosome Name"].append(i)
		i+=1
	#

	lineCount=0
	for line in infile:
		words=line.rstrip().split(",")
		# skip the chromosomes that are not 1..23 or X or Y
		ch1,ch2 = words[columnIndices["Chromosome Name"][0]],words[columnIndices["Chromosome Name"][1]]
		if ch1=="MT" or len(ch1)>2 or ch2=="MT" or len(ch2)>2:
			continue
		#
		gene1,gene2= words[columnIndices["GeneID1"]], words[columnIndices["GeneID2"]]
		homologyType= words[columnIndices["Homology Type"]]
		orthologyConfidence=words[columnIndices["Orthology confidence"]]
		percentIdentity=words[columnIndices["Percent Identity"]]
		# below checks ensure I only get one entry per g1-g2 pair
		if gene1 not in genePairsDic:
			genePairsDic[gene1]=[]
			genePairsDic[gene1].append([gene2,homologyType,orthologyConfidence,percentIdentity])
		else: 
			if gene2 not in [a[0] for a in genePairsDic[gene1]]:
				genePairsDic[gene1].append([gene2,homologyType,orthologyConfidence,percentIdentity])
			#
		#
		lineCount+=1
	#
	#print len(genePairsDic)
	types={}
	for g1 in genePairsDic:
		type=genePairsDic[g1][0][1] 
		num=len(genePairsDic[g1])
		if type=="ortholog_one2one" and num>1:
			sys.exit("Matching should be one2one but it isn't\t"+g1)
		if type not in types:
			types[type]=0
		types[type]+=1
	#
	infile.close()
	
	sys.stderr.write("Gene pair mappings summary: "+repr(types)+"\n\n")

	return genePairsDic


def parse_organism_GTF(orgID, infilename, outdir):
	"""
	This function parses a given Ensembl GTF file into the 
	internal data structure for that organism. 
	Does not output anything if outdir is "None"
	"""
	sys.stderr.write("Parsing organism GTF for "+orgID+" from file "+infilename+"\n")
	geneDic={}
	transcriptDic={}
	exonDic={}
	infoDic={} # build, version, accession
	if infilename.endswith(".gz"):
		infile=gzip.open(infilename,'rb')
	else:
		infile=open(infilename,'r')
	# parse the information header first
	elemCounts={"CDS" : 0, "exon" : 0, "gene" : 0, "start_codon" : 0, 
		"stop_codon" : 0, "transcript" : 0, "UTR" : 0, "other" : 0}
	lineCount=0
	lastReadExon="dummy"
	for line in infile:
		if line.startswith("#"):
			key, item = line.split()[:2]
			infoDic[key.split("-")[-1]]=item
		else:
			elemType=line.rstrip().split("\t")[2]
			chrName=line.rstrip().split("\t")[0]
			if chrName=="MT" or len(chrName)>2:
				#print chrName
				continue
			if elemType=="gene":
				newGene=EnsemblGene(line)
				geneDic[newGene.basicInfoDic["gene_id"]]=newGene
			elif elemType=="transcript":
				newTranscript=EnsemblTranscript(line)
				transcriptDic[newTranscript.basicInfoDic["transcript_id"]]=newTranscript
				geneId=newTranscript.basicInfoDic["gene_id"]
				geneDic[geneId].add_transcript(newTranscript)
			elif elemType=="exon":
				newExon=EnsemblExon(line)
				lastReadExon=newExon.basicInfoDic["exon_id"]
				
				# DELETEME!!
				#print("ALL_EXON_ENTRY\t%s\n" % (newExon.get_summary_string())),

				# Make sure to store certain information about the exon if it appears multiple times
				if lastReadExon not in exonDic:
					exonDic[lastReadExon]=newExon
				else:
					# DELETEME!!
					#print("DUPLICATE_EXON_ENTRY\t%s\t%s\n" % (exonDic[lastReadExon].get_summary_string(),newExon.get_summary_string())),
					exonDic[lastReadExon].add_another_instance(newExon)
				#
				geneId=newExon.basicInfoDic["gene_id"]
				transcriptId=newExon.basicInfoDic["transcript_id"]
				## Add the exon to a gene, don't care about ordering and simply owerwrite if exists
				geneDic[geneId].add_exon(newExon)
				## Add to exon to a transcript, make sure this exon insertion is ordered
				## Luckily entrys come ordered!
				## also same exon doesn't appear twice in once transcript so that case is not handled
				transcriptDic[transcriptId].add_exon(newExon)
				## no need for below line because Python is Pass-by-object-reference
				#geneDic[geneId].transcripts[transcriptId].add_exon(newExon)
			elif elemType=="CDS":
				#meaning previously read exon is completely/partially coding
				newLocus=EnsemblLocus(line)
				transcriptDic[transcriptId].handle_CDS(newLocus)
				exonDic[lastReadExon].handle_CDS(newLocus)
				# DELETEME!!
				#print "handleCDS\t%s\t%s\n" % (lastReadExon,newLocus.get_summary_string()),
			elif elemType=="UTR" or elemType=="stop_codon" or elemType=="start_codon":
				newLocus=EnsemblLocus(line)
				transcriptDic[transcriptId].add_locus(newLocus,elemType)
			#
			if elemType not in elemCounts:
				elemType="other"
			elemCounts[elemType]=elemCounts[elemType]+1
			lineCount+=1
			if lineCount%100000==0:
				sys.stderr.write(str(lineCount)+"\t")
		#
	#
	sys.stderr.write("\n")
	sys.stderr.write("GTF parsing summary: " +repr(elemCounts)+"\n\n")
	infile.close()
	if outdir!="None":
		print_some_summary(orgID, geneDic,transcriptDic,exonDic,elemCounts, outdir)
	return (geneDic,transcriptDic,exonDic,infoDic)


def print_some_summary(orgID, geneDic,transcriptDic,exonDic,elemCounts, outdir):
	"""
	Print a summary for the genes, transcripts and exons in the
	read GTF file.
	"""
	outfile=open(outdir+"/"+orgID+"-allGenes-GTFparsed.txt",'w')
	outfile.write("chrName\tstartCoord\tendCoord\tstrand\tgeneID\tgeneName\tgeneType\tnoOfTranscripts\tnoOfExons\telementType\n")
	for g in geneDic:
		outfile.write(geneDic[g].get_summary_string()+"\n")
		#print geneDic[g].get_summary_string()
	#
	outfile.close()

	totalNumberOfExons=0
	outfile=open(outdir+"/"+orgID+"-allTranscripts-GTFparsed.txt",'w')
	outfile.write("chrName\tstartCoord\tendCoord\tstrand\ttranscriptID\ttranscriptName\ttranscriptType\tgeneID\tgeneName\texonIDs\texonTypes\texonStarts\texonEnds\tcodingStarts\tcodingEnds\tstartCodon\tstopCodon\tUTRstarts\tUTRends\tproteinID\telementType\n")
	for t in transcriptDic:
		outfile.write(transcriptDic[t].get_summary_string()+"\n")
		totalNumberOfExons+=len(transcriptDic[t].exon_types)
		#print transcriptDic[t].get_summary_string()
	#
	outfile.close()
#	print ["totalNumberOfExons", totalNumberOfExons]

	outfile=open(outdir+"/"+orgID+"-allExons-GTFparsed.txt",'w')
	outfile.write("chrName\tstartCoord\tendCoord\tstrand\texonID\texonType\tcodingStart\tcodingEnd\ttranscriptIDs\texonNumbers\tgeneID\texonLength\tacceptor2bp\tdonor2bp\tavgCodingConsScore\tavgConsScore\tfirstMidLastCounts\telementType\n")
	for e in exonDic:
		outfile.write(exonDic[e].get_summary_string()+"\n")
		#print exonDic[e].get_summary_string()
	#
	outfile.close()

	return

def overlapping_combined( orig_data, reverse = False):
    """
    Return list of intervals with overlapping neighbours merged together
    Assumes sorted intervals unless reverse is set

    """
    if not orig_data or not len(orig_data): return []
    if len(orig_data) == 1:
        return orig_data

    new_data = []

    if reverse:
        data = orig_data[:]
        data.reverse()
    else:
        data = orig_data

    if not data[0][0] <= data[1][0]:
        print data, reverse
    assert(data[0][0] <= data[1][0])

    # start with the first interval
    prev_beg, prev_end = data[0]

    # check if any subsequent intervals overlap
    for beg, end in data[1:]:
        if beg - prev_end + 1 > 0:
            new_data.append((prev_beg, prev_end))
            prev_beg = beg
        prev_end = max(end, prev_end)

    new_data.append((prev_beg, prev_end))

    if reverse:
        new_data.reverse()
    return new_data


def get_overlap_between_intervals(a, b):
	"""
	Finds the overlap between two intervals end points inclusive.
	#Makes sure not to report overlap beyond either interval length.
	#	a=[10,20]; b=[10,20] --> f(a,b)=10  !(not 11)
	#	a=[10,20]; b=[20,30] --> f(a,b)=1
		a=[10,20]; b=[15,30] --> f(a,b)=6
	"""
	#lena=abs(float(a[1])-float(a[0]))
	#lenb=abs(float(b[1])-float(b[0]))
	overlap=max(0, min(float(a[1]), float(b[1])) - max(float(a[0]), float(b[0]))+1)
	#minlen=min(lena,lenb)
	#return min(minlen,overlap)
	return overlap

def sort_by_column(somelist, n):
	"""
	Given a list with 1 or more columns this functions sorts it according 
	to the desired column n [0 len(list)). Does this in-place.
	"""
	somelist[:] = [(x[n], x) for x in somelist]
	somelist.sort()
	somelist[:] = [val for (key, val) in somelist]
	return

def chr_name_conversion(chrIn,org):
	"""
	Given an identifier for the chromosome name (str) or a chromosome number (int) 
	and an organism this function converts the identifier to the other representation. 
	Example: 
		converts 'chrX' or 'X' to 23 for human
		converts 23 to 'chrX' 1 to 'chr1' for human
	"""
	if isinstance(chrIn, int): # int to str
		if org=='human':
			if	chrIn<23 and chrIn>0:
				chrOut='chr'+str(chrIn)
			elif chrIn==23:
				chrOut='chrX'
			elif chrIn==24:
				chrOut='chrY'
			else:
				return 'problem'
		elif org=='mouse':
			if	chrIn<20 and chrIn>0:
				chrOut='chr'+str(chrIn)
			elif chrIn==20:
				chrOut='chrX'
			elif chrIn==21:
				chrOut='chrY'
			else:
				return 'problem'
		else:
			chrOut='chr'+str(chrIn)
	else: # str to int
		if 'chr' in chrIn:
			chrIn=chrIn[:3] # cut the 'chr'
		if org=='human':
			if	chrIn=='X':
				chrOut=23
			elif chrIn=='Y':
				chrOut=24
			else:
				chrOut=int(chrIn)
		elif org=='mouse':
			if	chrIn=='X':
				chrOut=20
			elif chrIn=='Y':
				chrOut=21
			else:
				chrOut=int(chrIn)
	return chrOut



################################# BEGIN ExtendedExon ##################################
### NOT USED FOR NOW, NOT YET IMPLEMENTED ####
class ExtendedExon:
	"""
	This class is a container for exons that combines input from multiple different
	sources/files. Below is a list of these sources:
		- Ensembl Exon: This will initiate the instance of the ExtendedExon class 
		- LiftOver Files:
		- Genomedata Archive:
		- PhastCons Scores:
		- BLAT within species: 
	"""
	def __init__(self, ensemblExon):
		self.basicInfoDic= ensemblExon.basicInfoDic
################################# END ExtendedExon ##################################



################################# BEGIN EnsemblExon ##################################
class EnsemblExon:
	"""
	This class is a container for Ensembl exons
	"""
	def __init__(self, line):
		# parse the transcript line
		chr,d,elemType,start_coord,end_coord,d,strand,d,others=line.rstrip().split("\t")
		if elemType!="exon":
			sys.exit("Not an exon parsed in class EnsemblExon:\t"+elemType)
		#
		#basic information about the exon
		self.basicInfoDic={"chromosome" : chr, "start_coord" : int(start_coord), "end_coord" : int(end_coord), "strand" : strand}

		# there are 13 or more items for exons. We keep only 7 relevant ones.
		# e.g: gene_id "ENSG00000167468"; gene_version "14"; 
		#	transcript_id "ENST00000593032"; transcript_version "3"; exon_number "2"; 
		#	gene_name "GPX4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; 
		#	transcript_name "GPX4-006"; transcript_source "havana"; transcript_biotype "protein_coding"; 
		#	exon_id "ENSE00003420595"; exon_version "1"; tag "seleno"; tag "cds_end_NF"; 
		items=others.replace('"', '').split(";")
		for item in items:
			wds=item.lstrip().split()
			if len(wds)>1:
				key, val = wds[0],wds[1]
				if key in ["gene_id", "gene_name", "transcript_id", "transcript_name", "transcript_biotype", "exon_id", "exon_number"]:
					self.basicInfoDic[key]=val
		#

		self.exon_type = "nonCoding" # by default
		self.codingExon =  [-1,-1] # the coordinates of below codingExon will change if partialCoding or coding
		self.transcriptIds=[self.basicInfoDic["transcript_id"]]
		self.exonNumbers=[int(self.basicInfoDic["exon_number"])]
		self.acceptor2bp="NN"
		self.donor2bp="NN"
		self.phastConsScores=[]
		self.avgConsScore=0
		self.avgCodingConsScore=0
		self.firstMidLast=[0,0,0] # appearences of this exon as first, mid and last exons. Single exons are counted as first and last.
	#

	def handle_CDS(self,newLocus):
		lastSt, lastEn=self.basicInfoDic["start_coord"], self.basicInfoDic["end_coord"]
		newSt, newEn= newLocus.basicInfoDic["start_coord"], newLocus.basicInfoDic["end_coord"]	
		if lastSt==newSt and lastEn==newEn:
			self.exon_type="fullCoding"
		elif get_overlap_between_intervals([lastSt,lastEn], [newSt,newEn])>0:
			self.exon_type="partCoding"
		else:
			sys.exit("Reached a CDS entry that doesn't overlap with previous exon\t"\
				+newLocus.get_summary_string()+"\n")   
		#
		self.codingExon=[newSt,newEn]
	#
	def add_another_instance(self,newExon):
		# an exon may appear in only one gene but for many different transcripts
		self.transcriptIds.append(newExon.basicInfoDic["transcript_id"])
		self.exonNumbers.append(int(newExon.basicInfoDic["exon_number"]))
	#
	# data containers within this class
	__slots__ = ["basicInfoDic", "exon_type", "codingExon", "transcriptIds", "exonNumbers", \
			"acceptor2bp", "donor2bp", "phastConsScores", "avgCodingConsScore", "avgConsScore", "firstMidLast"]

	# get one liner summary of the given instance
	def get_summary_string(self):
		summary="chr"+self.basicInfoDic["chromosome"]+"\t"+str(self.basicInfoDic["start_coord"])+"\t"+\
			str(self.basicInfoDic["end_coord"])+"\t"+self.basicInfoDic["strand"]+"\t"+self.basicInfoDic["exon_id"]+"\t"+\
			self.exon_type +"\t"+str(self.codingExon[0])+"\t"+str(self.codingExon[1])+"\t"+\
			",".join(self.transcriptIds)+"\t"+",".join([str(e) for e in self.exonNumbers])+"\t"+\
			self.basicInfoDic["gene_id"]+"\t"+str(abs(self.basicInfoDic["end_coord"]-self.basicInfoDic["start_coord"]))+"\t"+\
			self.acceptor2bp+"\t"+self.donor2bp+"\t"+str(self.avgCodingConsScore)+"\t"+str(self.avgConsScore)+"\t"+\
			",".join([str(e) for e in self.firstMidLast])+"\texon"
			#self.basicInfoDic["transcript_id"]+"\t"+self.basicInfoDic["exon_number"]+"\t"+\
		return summary
################################# END EnsemblExon ##################################

################################# BEGIN EnsemblLocus ##################################
class EnsemblLocus:
	"""
	This class is a container for a basic locus that has chr, start, end, strand 
	fields. UTRs, start and stop codons from Ensembl are of this type.
	"""
	def __init__(self, line):
		# parse the locus line
		chr,d,elemType,start_coord,end_coord,d,strand,d,others=line.rstrip().split("\t")
		if elemType!="UTR" and elemType!="stop_codon" and elemType!="start_codon" and elemType!="CDS":
			sys.exit("Not a basic locus as intended parsed in from Ensemble line:\t"+line)
		#
		#basic information about the locus
		self.basicInfoDic={"chromosome" : chr, "start_coord" : int(start_coord), \
			"end_coord" : int(end_coord), "strand" : strand, "locus_type" : elemType}

		items=others.replace('"', '').split(";")
		for item in items:
			wds=item.lstrip().split()
			if len(wds)>1:
				key, val = wds[0],wds[1]
				if key in ["gene_id", "gene_name", "transcript_id", "transcript_name", "transcript_biotype", "protein_id"]:
					self.basicInfoDic[key]=val
		#


		#

	__slots__ = ["basicInfoDic"]

	# get one liner summary of the given instance
	def get_summary_string(self):
		summary="chr"+self.basicInfoDic["chromosome"]+"\t"+str(self.basicInfoDic["start_coord"])+"\t"+\
			str(self.basicInfoDic["end_coord"])+"\t"+self.basicInfoDic["strand"]+"\t"+self.basicInfoDic["locus_type"]+"\t"+\
			self.basicInfoDic["transcript_id"]+"\t"+self.basicInfoDic["transcript_name"]+"\t"+\
			self.basicInfoDic["transcript_biotype"]+"\t"+self.basicInfoDic["gene_id"]+"\t"+self.basicInfoDic["gene_name"]+"\tlocus"
		return summary
################################# END EnsemblLocus ##################################


################################# BEGIN EnsemblTranscript ##################################
class EnsemblTranscript:
	"""
	This class is a container for Ensembl transcripts
	"""
	def __init__(self, line):
		# parse the transcript line
		chr,d,elemType,start_coord,end_coord,d,strand,d,others=line.rstrip().split("\t")
		if elemType!="transcript":
			sys.exit("Not a transcript parsed in class EnsemblTranscript:\t"+elemType)
		#
		#basic information about the transcript	
		self.basicInfoDic={"chromosome" : chr, "start_coord" : int(start_coord), "end_coord" : int(end_coord), "strand" : strand}

		# there are 10 or more items for transcripts. We keep only 5 relevant ones.
		# e.g: gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; 
		#	transcript_version "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; 
		#	transcript_name "DDX11L1-002"; transcript_source "havana"; transcript_biotype "processed_transcript";
		#

		items=others.replace('"', '').split(";")
		for item in items:
			wds=item.lstrip().split()
			if len(wds)>1:
				key, val = wds[0],wds[1]
				#key, val = (item.lstrip().split())[0:2]
				if key in ["gene_id", "gene_name", "transcript_id", "transcript_name", "transcript_biotype"]:
					self.basicInfoDic[key]=val
		#
		if "gene_name" not in self.basicInfoDic: 
			self.basicInfoDic["gene_name"]=self.basicInfoDic["gene_id"]
		if "transcript_name" not in self.basicInfoDic: 
			self.basicInfoDic["transcript_name"]=self.basicInfoDic["transcript_id"]

		self.start_codon=["cds_start_NF"] # by default make them non-confirmed
		self.stop_codon=["cds_stop_NF"]  # by default make them non-confirmed
		self.exons = []
		self.codingExons = []
		self.exon_types= []
		self.protein_id="None"
		self.UTRs=[]
	#

	# adding an exon to the list of exons of the transcript "in order"
	def add_exon(self,newExon):
		exonCountSoFar=len(self.exons)
		#print "to add\t"+str(newExon.basicInfoDic["exon_number"])+"\t"+str(newExon.basicInfoDic["exon_id"])
		if int(newExon.basicInfoDic["exon_number"])==exonCountSoFar+1:
			exonEntry=[newExon.basicInfoDic["start_coord"],newExon.basicInfoDic["end_coord"],
				newExon.basicInfoDic["exon_id"]]
			self.exons.append(exonEntry)
 			# the coordinates of below codingExon will change if partialCoding or coding
			self.codingExons.append([-1,-1])
			exonType="nonCoding" # by default
			self.exon_types.append(exonType)
		else:
			sys.exit("Exon entry is being entered out of order to the transcript\t"
				+self.basicInfoDic["transcript_id"])
		#
	#
	def add_locus(self,newLocus,locus_type):
		if locus_type=='start_codon':
			self.start_codon=[newLocus.basicInfoDic["start_coord"],newLocus.basicInfoDic["end_coord"]]
		elif locus_type=='stop_codon':
			self.stop_codon=[newLocus.basicInfoDic["start_coord"],newLocus.basicInfoDic["end_coord"]]
		elif locus_type=='UTR':
			self.UTRs.append([newLocus.basicInfoDic["start_coord"],newLocus.basicInfoDic["end_coord"]])
		else:
			sys.exit("Unknow locus type being inserted to transcript\t"\
				+self.basicInfoDic["transcript_id"])   
	#

	def handle_CDS(self,newLocus):
		exonType="nonCoding" # by default
		exonCountSoFar=len(self.exons)
		lastAddedExon=self.exons[exonCountSoFar-1]
		lastSt,lastEn=self.exons[exonCountSoFar-1][0:2]
		newSt, newEn= newLocus.basicInfoDic["start_coord"], newLocus.basicInfoDic["end_coord"]		
		if lastSt==newSt and lastEn==newEn:
			exonType="fullCoding"
		elif get_overlap_between_intervals([lastSt,lastEn], [newSt,newEn])>0:
			exonType="partCoding"
		else:
			sys.exit("Reached a CDS entry that doesn't overlap with previous exon\t"\
				+newLocus.get_summary_string()+"\n")
		#
		self.codingExons[exonCountSoFar-1]=[newSt,newEn]
		self.exon_types[exonCountSoFar-1]=exonType # replace with the previous nonCoding tag
		self.protein_id=newLocus.basicInfoDic["protein_id"]
	#

	# data containers within this class
	__slots__ = [
			"basicInfoDic", 
			"start_codon",
			"stop_codon",
			"exons",
			"codingExons",
			"exon_types",
			"protein_id",
			"UTRs"
	]
	# get one liner summary of the given instance
	def get_summary_string(self):
		if len(self.UTRs)==0:
			self.UTRs.append(["None","None"])
		#
		summary="chr"+self.basicInfoDic["chromosome"]+"\t"+str(self.basicInfoDic["start_coord"])+"\t"+\
			str(self.basicInfoDic["end_coord"])+"\t"+self.basicInfoDic["strand"]+"\t"+self.basicInfoDic["transcript_id"]+"\t"+\
			self.basicInfoDic["transcript_name"]+"\t"+self.basicInfoDic["transcript_biotype"]+"\t"+\
			self.basicInfoDic["gene_id"]+"\t"+self.basicInfoDic["gene_name"]+"\t"+\
			",".join([e[2] for e in self.exons])+"\t"+",".join(self.exon_types)+"\t"+\
			",".join([str(e[0]) for e in self.exons])+"\t"+",".join([str(e[1]) for e in self.exons])+"\t"+\
			",".join([str(e[0]) for e in self.codingExons])+"\t"+",".join([str(e[1]) for e in self.codingExons])+"\t"+\
			",".join([str(i) for i in self.start_codon])+"\t"+",".join([str(i) for i in self.stop_codon])+"\t"+\
			",".join([str(e[0]) for e in self.UTRs])+"\t"+",".join([str(e[1]) for e in self.UTRs])+"\t"+\
			self.protein_id+"\t"+"transcript"
		return summary
################################# END EnsemblTranscript ##################################


			
################################# BEGIN EnsemblGene ##################################
class EnsemblGene:
	"""
	This class is a container for Ensembl genes.
	"""
	def __init__(self, line):
		# parse the gene line
		chr,d,elemType,start_coord,end_coord,d,strand,d,others=line.rstrip().split("\t")
		if elemType!="gene":
			sys.exit("Not a gene parsed in class EnsemblGene:\t"+elemType)
		#
		#basic information about the gene	
		self.basicInfoDic={"chromosome" : chr, "start_coord" : int(start_coord), "end_coord" : int(end_coord), "strand" : strand}

		# there are 5 items for genes: 
		# e.g: gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed";
		items=others.replace('"', '').split(";")
		for item in items:
			if len(item)>1:
				key, val = item.lstrip().split()
				self.basicInfoDic[key]=val
		if "gene_name" not in self.basicInfoDic: 
			self.basicInfoDic["gene_name"]=self.basicInfoDic["gene_id"]
		self.exons = {}
		self.transcripts = {}
	#

	# data containers within this class
	__slots__ = [
			"basicInfoDic", 
			"exons",
			"transcripts"]

	# adding a transcript to the gene
	def add_transcript(self,newTranscript):
		self.transcripts[newTranscript.basicInfoDic["transcript_id"]]=newTranscript
	# adding an exon to the list of exons of the gene
	def add_exon(self,newExon):
		self.exons[newExon.basicInfoDic["exon_id"]]=\
			[newExon.basicInfoDic["start_coord"],newExon.basicInfoDic["end_coord"]]
	# get one liner summary of the given instance
	def get_summary_string(self):
		summary="chr"+self.basicInfoDic["chromosome"]+"\t"+str(self.basicInfoDic["start_coord"])+"\t"+\
			str(self.basicInfoDic["end_coord"])+"\t"+self.basicInfoDic["strand"]+"\t"+self.basicInfoDic["gene_id"]+"\t"+\
			self.basicInfoDic["gene_name"]+"\t"+self.basicInfoDic["gene_biotype"]+"\t"+\
			str(len(self.transcripts))+"\t"+str(len(self.exons))+"\tgene"
		return summary
	#
	#def gene_wrap_up():
	#	self.beg = min(self.exons[e][0] for e in self.exons)
	#	self.end = max(self.exons[e][1] for e in self.exons)
	#

################################# END EnsemblGene ##################################


def convert_UCSC_to_bed_format(l):
	"""
	Given a locus in UCSC format this function converts it to bed format with 3 fields
	chr1:121-21111 --> ['chr1', 121, 21111]
	"""
	chr=l[:l.find(':')]
	st=int(l[l.find(':')+1:l.find('-')])
	en=int(l[l.find('-')+1:])
	return (chr,st,en)


def consistency_check(org1TOorg2,org2TOorg1):
	"""
	Check the consistency between the two matchings (e.g. human-to-mouse, mouse-to-human)
	read from separate Ensembl file. This function will do nothing if all is consistent.
	"""
	for g1 in org1TOorg2:
		type=org1TOorg2[g1][0][1]
		if type=="ortholog_one2one":
			g2=org1TOorg2[g1][0][0]
			if g2 not in org2TOorg1:
				sys.exit("Reverse entry for a one2one match couldn't be found\t"+g1+"\t"+g2)
			elif org2TOorg1[g2][0][0]!=g1:
				sys.exit("Reverse entry for a one2one match mismatches with original one\t"+g1+"\t"+g2)
			# else good
		else:
			for oneMatch1 in org1TOorg2[g1]:
				g2=oneMatch1[0]
				if g2 not in org2TOorg1:
					sys.exit("Reverse entry for a NON-one2one match couldn't be found\t"+g1+"\t"+g2)
				else:
					reverseFound=False
					for oneMatch2 in org2TOorg1[g2]:
						if oneMatch2[0]==g1:
							reverseFound=True
							break
					#
					if reverseFound==False:
						sys.exit("Reverse entry for a NON-one2one match mismatches with original one\t"+g1+"\t"+g2)
					# else good
				#
	#
	for g1 in org2TOorg1:
		type=org2TOorg1[g1][0][1]
		if type=="ortholog_one2one":
			g2=org2TOorg1[g1][0][0]
			if g2 not in org1TOorg2:
				sys.exit("Reverse entry for a one2one match couldn't be found\t"+g1+"\t"+g2)
			elif org1TOorg2[g2][0][0]!=g1:
				sys.exit("Reverse entry for a one2one match mismatches with original one\t"+g1+"\t"+g2)
			# else good
		else:
			for oneMatch1 in org2TOorg1[g1]:
				g2=oneMatch1[0]
				if g2 not in org1TOorg2:
					sys.exit("Reverse entry for a NON-one2one match couldn't be found\t"+g1+"\t"+g2)
				else:
					reverseFound=False
					for oneMatch2 in org1TOorg2[g2]:
						if oneMatch2[0]==g1:
							reverseFound=True
							break
					#
					if reverseFound==False:
						sys.exit("Reverse entry for a NON-one2one match mismatches with original one\t"+g1+"\t"+g2)
					# else good
				#
	#
	return


def pickle_one2one_genePairs_allInfo(genePairsDic,geneDic1,geneDic2,exonDic1,exonDic2,transcriptDic1,transcriptDic2,outdir):
	"""
	Pickle the gene, transcript and exon dictionaries for each pair of ortholog_one2one genes.
	There are around 16.5k such genes for human-mouse and 15.8k are protein_coding pairs.
	"""
	geneOnlyOrthoDic1,transcriptOnlyOrthoDic1,exonOnlyOrthoDic1={},{},{}
	geneOnlyOrthoDic2,transcriptOnlyOrthoDic2,exonOnlyOrthoDic2={},{},{}

	for g1 in genePairsDic:
		type=genePairsDic[g1][0][1]
		if type=="ortholog_one2one":
			g2=genePairsDic[g1][0][0]
		else:
			continue

		if geneDic1[g1].basicInfoDic["gene_biotype"]!="protein_coding" or geneDic2[g2].basicInfoDic["gene_biotype"]!="protein_coding":
			continue

		# small dictionaries that have only the relevant stuff for one gene pair
		newGeneDic1={}; newGeneDic2={}
		newExonDic1={}; newExonDic2={}
		newTranscriptDic1={}; newTranscriptDic2={}
		#
		newGeneDic1[g1]=geneDic1[g1]
		newGeneDic2[g2]=geneDic2[g2]
		geneOnlyOrthoDic1[g1]=geneDic1[g1]
		geneOnlyOrthoDic2[g2]=geneDic2[g2]
		for tId in geneDic1[g1].transcripts:
			newTranscriptDic1[tId]=transcriptDic1[tId]
			transcriptOnlyOrthoDic1[tId]=transcriptDic1[tId]
		for tId in geneDic2[g2].transcripts:
			newTranscriptDic2[tId]=transcriptDic2[tId]
			transcriptOnlyOrthoDic2[tId]=transcriptDic2[tId]
		#
		for eId in geneDic1[g1].exons:
			newExonDic1[eId]=exonDic1[eId]
			exonOnlyOrthoDic1[eId]=exonDic1[eId]
		for eId in geneDic2[g2].exons:
			newExonDic2[eId]=exonDic2[eId]
			exonOnlyOrthoDic2[eId]=exonDic2[eId]
		#
		os.system("mkdir -p "+ outdir+"/"+g1+"-"+g2)
		#print geneDic1[g1].get_summary_string()+"\t"+geneDic2[g2].get_summary_string()
		outfilename=outdir+"/"+g1+"-"+g2+"/org1.pickledDictionaries"
		pickle.dump((newGeneDic1,newTranscriptDic1,newExonDic1), open(outfilename,"wb"))
		outfilename=outdir+"/"+g1+"-"+g2+"/org2.pickledDictionaries"
		pickle.dump((newGeneDic2,newTranscriptDic2,newExonDic2), open(outfilename,"wb"))
		# to load use:
		#geneDic1,transcriptDic1,exonDic1=pickle.load(open("pickled.stuff","rb"))
	#
	
	return (geneOnlyOrthoDic1,transcriptOnlyOrthoDic1,exonOnlyOrthoDic1, geneOnlyOrthoDic2,transcriptOnlyOrthoDic2,exonOnlyOrthoDic2)

def print_one2one_genePairs(genePairsDic, geneDic1,geneDic2,outfilename):
	"""
	Print one liner for each pair of genes that match each other one to one.
	There are around 16.5k such genes for human-mouse and 15.8k are protein_coding pairs.
	"""
	outfile=open(outfilename,'w')
	outfile.write("chrName1\tstart_coord1\tend_coord1\tstrand1\tgeneID1\tgeneName1\tgeneType1\tnoOfTranscripts1\tnoOfExons1\ttype1\t")
	outfile.write("chrName2\tstart_coord2\tend_coord2\tstrand2\tgeneID2\tgeneName2\tgeneType2\tnoOfTranscripts2\tnoOfExons2\ttype2\n")
	#print "chrName1\tstart_coord1\tend_coord1\tstrand1\tgeneID1\tgeneName1\tgeneType1\tnoOfTranscripts1\tnoOfExons1\ttype1\t",
	#print "chrName2\tstart_coord2\tend_coord2\tstrand2\tgeneID2\tgeneName2\tgeneType2\tnoOfTranscripts2\tnoOfExons2\ttype2"
	for g1 in genePairsDic:
		type=genePairsDic[g1][0][1]
		if type=="ortholog_one2one":
			g2=genePairsDic[g1][0][0]
		else:
			continue
		outfile.write(geneDic1[g1].get_summary_string()+"\t"+geneDic2[g2].get_summary_string()+"\n")
		#print geneDic1[g1].get_summary_string()+"\t"+geneDic2[g2].get_summary_string()
	#
	outfile.close()
	return

def print_one2one_transcriptListPairs(genePairsDic, geneDic1,geneDic2,transcriptDic1,transcriptDic2,orgId1,orgId2,outdir):
	"""
	Print the lists of transcripts for each one to one mapped gene pair. 
	There are around 16.5k such genes for human-mouse and 15.8k are protein_coding pairs.
	"""
	for g1 in genePairsDic:
		type=genePairsDic[g1][0][1]
		if type=="ortholog_one2one":
			g2=genePairsDic[g1][0][0]
		else:
			continue
		#
		outdirTemp=outdir+"/"+g1+"-"+g2; os.system("mkdir -p "+outdirTemp)
		outfile1=open(outdirTemp+"/"+orgId1+"_transcripts.bed",'w')
		outfile2=open(outdirTemp+"/"+orgId2+"_transcripts.bed",'w')
		
		transcripts1=geneDic1[g1].transcripts
		transcripts2=geneDic2[g2].transcripts
		for t1 in transcripts1:
			outfile1.write(transcriptDic1[t1].get_summary_string()+"\n")
		for t2 in transcripts2:
			outfile2.write(transcriptDic2[t2].get_summary_string()+"\n")
		#
		outfile1.close()
		outfile2.close()
	#
	return

def print_one2one_exonListPairs(genePairsDic, geneDic1,geneDic2,exonDic1,exonDic2,orgId1,orgId2,outdir):
	"""
	Print the lists of exons for each one to one mapped gene pair. 
	There are around 16.5k such genes for human-mouse and 15.8k are protein_coding pairs.
	"""
	for g1 in genePairsDic:
		type=genePairsDic[g1][0][1]
		if type=="ortholog_one2one":
			g2=genePairsDic[g1][0][0]
		else:
			continue
		#
		outdirTemp=outdir+"/"+g1+"-"+g2; os.system("mkdir -p "+outdirTemp)
		outfile1=open(outdirTemp+"/"+orgId1+"_exons.bed",'w')
		outfile2=open(outdirTemp+"/"+orgId2+"_exons.bed",'w')
		
		exons1=geneDic1[g1].exons
		exons2=geneDic2[g2].exons
		for e1 in exons1:
			outfile1.write(exonDic1[e1].get_summary_string()+"\n")
		for e2 in exons2:
			outfile2.write(exonDic2[e2].get_summary_string()+"\n")
		#
	#	print exonDic[e].get_summary_string()

		outfile1.close()
		outfile2.close()
		
	#
	return

def extract_fasta_files_for_exons(refGD,exonDic,typ,fivePrimeFlank,threePrimeFlank,outfilename):
	"""
	With the help of genomedata archive extract the nucleotide sequences 
	from and around each exon and write them in a .fa file.
	refGD is the genomedata archive created for the reference genome.
	typ can be one of the following:
		"allExon": Extract the sequence of the whole exon.
		"allExonPlusMinus": Like allExon but with flanking 5' and 3'.
		"intronExon": Extract the sequence from the juction of this
		   exon and the previous intron.
		"exonIntron": Extract the sequence from the juction of this
		   exon and the next intron.
	fivePrimeFlank is the amount to extract extra from the 5' end.
	threePrimeFlank is the amount to extract extra from the 3' end.

	"""
	if typ=="allExon":
		fivePrimeFlank=0; threePrimeFlank=0
	#
	outfile=open(outfilename,'w')
	with Genome(refGD) as genome:
		for id in exonDic:
			e=exonDic[id]
			ch,st,en="chr"+e.basicInfoDic["chromosome"], e.basicInfoDic["start_coord"], e.basicInfoDic["end_coord"]
			strand,id=e.basicInfoDic["strand"], e.basicInfoDic["exon_id"]
			# off by one error fix by -1
			st=int(st)-1
			en=int(en)-1
			if strand=="+":
				if typ=="intronExon":
					en=st # make sure we're around the first bp of exon
					st=st-fivePrimeFlank # make sure 5' part is of size fivePrimeFlank including st
					en=en+threePrimeFlank # make sure 3' part is of size threePrimeFlank including en
				elif typ=="exonIntron":
					st=en # make sure we're around the last bp of exon
					st=st-fivePrimeFlank+1
					en=en+threePrimeFlank+1
				elif typ=="allExonPlusMinus" or typ=="allExon":
					st=st-fivePrimeFlank
					en=en+threePrimeFlank+1
				#
				id=id+"_plusStrand"
				sq=genome[ch].seq[st:en].tostring().lower().upper()
			else:
				if typ=="intronExon":
					st=en # make sure we're around the first bp of exon
					en=en+fivePrimeFlank+1
					st=st-threePrimeFlank+1
				elif typ=="exonIntron":
					en=st # make sure we're around the last bp of exon
					en=en+fivePrimeFlank
					st=st-threePrimeFlank
				elif typ=="allExonPlusMinus" or typ=="allExon":
					st=st-threePrimeFlank
					en=en+fivePrimeFlank+1
				#
				id=id+"_minusStrand"
				sq=genome[ch].seq[st:en].tostring()
				#sq=sq.lower()[::-1].upper() # reverse
				#sq=sq.lower().translate(complement).upper() # complement
				sq=sq.lower().translate(complement)[::-1].upper() # reverse complement
			#
			outfile.write(">"+id+"_"+typ+"\n")
			outfile.write(sq+"\n")
		#
	#
	outfile.close()
	return

def extract_conservation_stats_for_exons(refGD,exonDic,typ,fivePrimeFlank,threePrimeFlank,outfilename):
	"""
	With the help of genomedata archive extract the nucleotide sequences 
	from and around each exon and convservation scores and write them to a file.
	refGD is the genomedata archive created for the reference genome.
	typ can be one of the following:
		"allExon": Extract the sequence of the whole exon.
		"allExonPlusMinus": Like allExon but with flanking 5' and 3'.
		"intronExon": Extract the sequence from the juction of this
		   exon and the previous intron.
		"exonIntron": Extract the sequence from the juction of this
		   exon and the next intron.
	fivePrimeFlank is the amount to extract extra from the 5' end.
	threePrimeFlank is the amount to extract extra from the 3' end.
	IF outfilename is "None" then no output file is written, only 
	relevant fields are added to the exon in exonDic.
	"""
	sys.stderr.write("Extracting conservation stats and acceptor donor sites for exons from genomedata archive\n")
	# this is the trackname for phastCons scores loaded from wig files
	trackName="phastCons"
	#

	if typ=="allExon":
		fivePrimeFlank=0; threePrimeFlank=0
	#
	if outfilename!="None":
		outfile=open(outfilename,'w')
		# header line
		outfile.write("CHR\tstart\tend\tstrand\tExonID\tacceptor2bp\tdonor2bp\tpreAcceptorCons\taccepterCons1\taccepterCons2\texon5primeCons\texonMidCons\texon3primeCons\tdonorCons1\tdonorCons2\tpostDonorCons\n")
	#
	with Genome(refGD) as genome:
		lineCount=0
		for id in exonDic:
			e=exonDic[id]
			ch,st,en="chr"+e.basicInfoDic["chromosome"], e.basicInfoDic["start_coord"], e.basicInfoDic["end_coord"]
			if e.exon_type=="partCoding":
				codingSt=min(int(e.codingExon[0])-1,int(e.codingExon[1])-1)
				codingEn=max(int(e.codingExon[0])-1,int(e.codingExon[1])-1)
			#
			stOrig=st; enOrig=en;
			strand,id=e.basicInfoDic["strand"], e.basicInfoDic["exon_id"]
			# off by one error fix by -1
			st=int(st)-1
			en=int(en)-1
			if strand=="+":
				if typ=="intronExon":
					en=st # make sure we're around the first bp of exon
					st=st-fivePrimeFlank # make sure 5' part is of size fivePrimeFlank including st
					en=en+threePrimeFlank # make sure 3' part is of size threePrimeFlank including en
				elif typ=="exonIntron":
					st=en # make sure we're around the last bp of exon
					st=st-fivePrimeFlank+1
					en=en+threePrimeFlank+1
				elif typ=="allExonPlusMinus":
					st=st-fivePrimeFlank
					en=en+threePrimeFlank+1
				#
				#id=id+"_plusStrand"
				sq=genome[ch].seq[st:en].tostring().lower().upper()
				allScores=(genome[ch])[st:en,trackName]
			else:
				if typ=="intronExon":
					st=en # make sure we're around the first bp of exon
					en=en+fivePrimeFlank+1
					st=st-threePrimeFlank+1
				elif typ=="exonIntron":
					en=st # make sure we're around the last bp of exon
					en=en+fivePrimeFlank
					st=st-threePrimeFlank
				elif typ=="allExonPlusMinus":
					st=st-threePrimeFlank
					en=en+fivePrimeFlank+1
				#
				#id=id+"_minusStrand"
				sq=genome[ch].seq[st:en].tostring()
				sq=sq.lower().translate(complement)[::-1].upper() # reverse complement
				allScores=(genome[ch])[st:en,trackName][::-1]
			#
			if e.exon_type=="partCoding":
				codingScores=(genome[ch])[codingSt:codingEn,trackName][::-1]
			### Extract all the scores to be written to the output file ###
			acceptor2bp=sq[fivePrimeFlank-2:fivePrimeFlank]
			donor2bp=(sq[-threePrimeFlank:])[0:2]
			#
			x=allScores[:fivePrimeFlank-2]
			preAcceptorCons=np.nanmean(x)
			#
			accepterCons1=allScores[fivePrimeFlank-2]
			accepterCons2=allScores[fivePrimeFlank-1]
			#
			x=allScores[fivePrimeFlank:fivePrimeFlank+(fivePrimeFlank-2)]
			exon5primeCons=np.nanmean(x)
			#
			x=allScores[fivePrimeFlank+(fivePrimeFlank-2):-(threePrimeFlank+(threePrimeFlank-2))]
			exonMidCons=np.nanmean(x)
			#
			x=allScores[-(threePrimeFlank+(threePrimeFlank+2)):-threePrimeFlank]
			exon3primeCons=np.nanmean(x)
			#
			donorCons1=allScores[-threePrimeFlank]
			donorCons2=allScores[-threePrimeFlank+1]
			#
			x=allScores[-threePrimeFlank+2:]
			postDonorCons=np.nanmean(x)
			#
			#first20bp=allScores[:20]
			#outfile.write("%s\t%d\t%d\t%s\t%s\t%s\t%s\t" % (ch,stOrig,enOrig,strand,id,acceptor2bp,donor2bp))
			#outfile.write("\t".join([repr(x) for x in first20bp])+"\n")
			exonDic[id].acceptor2bp=acceptor2bp
			exonDic[id].donor2bp=donor2bp
			exonDic[id].phastConsScores=allScores
			exonDic[id].avgConsScore=np.nanmean(allScores[fivePrimeFlank:-threePrimeFlank])
			if e.exon_type=="partCoding":
				exonDic[id].avgCodingConsScore=np.nanmean(codingScores)

			if lineCount%100000==0:
				sys.stderr.write(str(lineCount)+"\t")
			lineCount+=1

			if outfilename!="None":
				outfile.write("%s\t%d\t%d\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n" % \
					(ch,stOrig,enOrig,strand,id,acceptor2bp,donor2bp, preAcceptorCons, accepterCons1, accepterCons2, exon5primeCons,\
					 exonMidCons, exon3primeCons, donorCons1, donorCons2, postDonorCons, exonDic[id].avgCodingConsScore, exonDic[id].avgConsScore) )
			#
			###
		#
	#
	sys.stderr.write("\n\n")
	if outfilename!="None":
		outfile.close()
	return exonDic

def assign_firstMidLast_exon_counts(exonDic,transcriptDic):
	for e in exonDic:
		for i in range(len(exonDic[e].exonNumbers)):
			transcriptLength=len(transcriptDic[exonDic[e].transcriptIds[i]].exons)
			tempExNo=exonDic[e].exonNumbers[i]
			#print [transcriptLength, tempExNo]
			#single exon
			if tempExNo==1 and tempExNo==transcriptLength:
				exonDic[e].firstMidLast[0]+=1
				exonDic[e].firstMidLast[2]+=1
			#first exon
			elif tempExNo==1:
				exonDic[e].firstMidLast[0]+=1
			#last exon
			elif tempExNo==transcriptLength:
				exonDic[e].firstMidLast[2]+=1
			else:
				exonDic[e].firstMidLast[1]+=1
		#
		#if exonDic[e].firstMidLast[0]>0 or exonDic[e].firstMidLast[2]>0:
		#print exonDic[e].get_summary_string()
	#
	return exonDic

#   Testing functionalities
def main(argv):
	orgId1="human"; orgId2="mouse";
	refGD1="/home/fao150/proj/2015orthoR01/results/2015-03-17_creating-genomedata-archives-for-refs/hg38"
	refGD2="/home/fao150/proj/2015orthoR01/results/2015-03-17_creating-genomedata-archives-for-refs/mm10"

#	outdir="GTFsummaries"; 
	if len(argv)==1:
		return

	outdir=argv[1]
	os.system("mkdir -p "+outdir)
		
	infilename="/projects/b1017/shared/Ensembl-files/Homo_sapiens.GRCh38.78.gtf.gz"
	geneDic1,transcriptDic1,exonDic1,infoDic1=parse_organism_GTF(orgId1, infilename, outdir)
	
	infilename="/projects/b1017/shared/Ensembl-files/Mus_musculus.GRCm38.78.gtf.gz"
	geneDic2,transcriptDic2,exonDic2,infoDic2=parse_organism_GTF(orgId2, infilename, outdir)

	## these two files were downloaded by hand selecting columns from Ensembl's Biomart
	## I weren't able to redo the same column selections recently so I decided to switch to
	## parsing the orthology information from readily available Ensembl files like below ones:
	## ftp://ftp.ensembl.org/pub/release-80/mysql/ensembl_mart_80/
	## hsapiens_gene_ensembl__homolog_mmus__dm.txt.gz
	#infilename="/projects/b1017/shared/Ensembl-files/Ensembl-human-GRCh38-to-mouse-GRCm38.p3.txt.gz"
	#genePairsHumanToMouse=parse_ensembl_gene_pairings(orgId1,orgId2,infilename)
	#infilename="/projects/b1017/shared/Ensembl-files/Ensembl-mouse-GRCm38.p3-to-human-GRCh38.txt.gz"
	#genePairsMouseToHuman=parse_ensembl_gene_pairings(orgId2,orgId1,infilename)
	#consistency_check(genePairsHumanToMouse,genePairsMouseToHuman)
	## if consistency check is ok then just use one side. This is OK for one2one mappings.
	#genePairsDic=genePairsHumanToMouse
	#pickle_one2one_genePairs_allInfo(genePairsDic,geneDic1,geneDic2,exonDic1,exonDic2,transcriptDic1,transcriptDic2,outdir) 
	
	infilename="/projects/b1017/shared/Ensembl-files/hsapiens_gene_ensembl__homolog_mmus__dm.txt.gz"
	proteinToGeneDic,genePairsDic,proteinPairsDic=parse_ensembl_geneAndProtein_pairings(infilename,{},{})
	print ["1",len(proteinToGeneDic),len(genePairsDic),len(proteinPairsDic)]

	infilename="/projects/b1017/shared/Ensembl-files/mmusculus_gene_ensembl__homolog_hsap__dm.txt.gz"
	proteinToGeneDic,genePairsDic,proteinPairsDic=parse_ensembl_geneAndProtein_pairings(infilename,proteinToGeneDic,proteinPairsDic)
	print ["2",len(proteinToGeneDic),len(genePairsDic),len(proteinPairsDic)]

	
	exonDic1=assign_firstMidLast_exon_counts(exonDic1,transcriptDic1)
	exonDic2=assign_firstMidLast_exon_counts(exonDic2,transcriptDic2)

	typ="allExonPlusMinus"
	outfilename="None"
	fivePrimeFlank=12; threePrimeFlank=12
	exonDic1=extract_conservation_stats_for_exons(refGD1,exonDic1,typ,fivePrimeFlank,threePrimeFlank,outfilename)
	exonDic2=extract_conservation_stats_for_exons(refGD2,exonDic2,typ,fivePrimeFlank,threePrimeFlank,outfilename)

	outdir=argv[1]+"/after"
	os.system("mkdir -p "+outdir)
	print_some_summary(orgId1, geneDic1,transcriptDic1,exonDic1,{}, outdir)
	print_some_summary(orgId2, geneDic2,transcriptDic2,exonDic2,{}, outdir)

#	outdir="perGenePairExonLists"
	if len(argv)==2:
		return

	outdir=argv[2]
	os.system("mkdir -p "+outdir)

	pickle_one2one_genePairs_allInfo(genePairsDic,geneDic1,geneDic2,exonDic1,exonDic2,transcriptDic1,transcriptDic2,outdir) 

#	outfilename=outdir+"/genePairsSummary-one2one.txt"
#	print_one2one_genePairs(genePairsDic,geneDic1,geneDic2,outfilename) # either way is ok since one2one
	#
#	print_one2one_exonListPairs(genePairsDic,geneDic1,geneDic2,exonDic1,exonDic2,orgId1,orgId2,outdir)
#	print_one2one_transcriptListPairs(genePairsDic,geneDic1,geneDic2,transcriptDic1,transcriptDic2,orgId1,orgId2,outdir)

	return

if __name__ == "__main__":
	main(sys.argv)


