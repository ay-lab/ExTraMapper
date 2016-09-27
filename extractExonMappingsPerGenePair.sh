#!/bin/bash -ex
set -o pipefail
set -o errexit 

### This script extracts the exon mappings for a given gene pair from all three methods and outputs them ###
# Method1: Zhang-GenomeBiology-2009 
# Method2: Fu-BMCGemomics-2012 (Orthoexon) 
# Method3: ExTraMapper (this paper)

# read all the configuration parameters first
source config.conf
export EXTRAMAPPER_DIR=$EXTRAMAPPERDIR

#gene1=ENSG00000157764; gene2=ENSMUSG00000002413; # BRAF-Braf
gene1=$1 # user input
gene2=$2 # user input

dataDir=${EXTRAMAPPER_DIR}/preprocess/data
dataDirPerPair=${EXTRAMAPPER_DIR}/preprocess/data/$org1-$org2
referenceGenomesDir=$dataDir/reference_genomes
chainsDir=$dataDir/liftover_chains
ensemblDir=$dataDirPerPair/ensemblDownloads
genomedataDir=$dataDirPerPair/genomedataArchives
phastConsDir=$dataDirPerPair/phastConsScores
GTFsummaryDir=$dataDirPerPair/GTFsummaries
perGenePairPickleDir=$dataDirPerPair/perGenePairPickledInfo
liftOverFilesDir=$dataDirPerPair/liftoverRelatedFiles
perExonLiftoverDir=$dataDirPerPair/perExonLiftoverCoords
#


cat $GTFsummaryDir/org1-allExons-GTFparsed.txt | grep $gene1 |  awk '{print $5"\t"$0}' | sort -k1,1 > myg1.tmp
cat $GTFsummaryDir/org2-allExons-GTFparsed.txt | grep $gene2 |  awk '{print $5"\t"$0}' | sort -k1,1 > myg2.tmp

for paper in Zhang-GenomeBiology-2009 Fu-BMCGemomics-2012 ExTraMapper; do
#	paper=Fu-BMCGemomics-2012 # either Fu-BMCGemomics-2012 or Zhang-GenomeBiology-2009
	orthoExonFile=$paper/Human_and_mouse_orthologous_exons-Ensembl-hg38-to-mm10.txt

	if [[ $paper == "ExTraMapper" ]]; then
		orthoExonFile=output/human-mouse/$gene1-$gene2/exonLevelMappings-1.0.txt
	fi

	cat $orthoExonFile | awk '{print $5"\t"$0}' | sort -k1,1 > oeg1.tmp
	cat $orthoExonFile | awk '{print $11"\t"$0}' | sort -k1,1 > oeg2.tmp
	cat $orthoExonFile | awk '{print $5"-"$11}' | sort > oepairs.tmp

	join oeg1.tmp myg1.tmp > int1.tmp
	join oeg2.tmp myg2.tmp > int2.tmp

	cat int1.tmp | awk '{print $6"-"$12,$19,$22,$24,$25,$26,$27,$28,$29,$30}' | sort -k1,1 > int1.sorted.tmp
	cat int2.tmp | awk '{print $6"-"$12,$19,$22,$24,$25,$26,$27,$28,$29,$30}' | sort -k1,1 > int2.sorted.tmp
	echo -ne "exonPair\texonType1\ttranscriptIDs1\tgeneID1\texonLength1\tacceptor2bp1\tdonor2bp1\tavgCodingConsScore1\tavgConsScore1\tfirstMidLastCounts1\texonType2\ttranscriptIDs2\tgeneID2\texonLength2\tacceptor2bp2\tdonor2bp2\tavgCodingConsScore2\tavgConsScore2\tfirstMidLastCounts2\n" > intint.tmp
	join int1.sorted.tmp int2.sorted.tmp | sed 's/ /\t/g' | sort -k1,1 >> intint.tmp

	mv intint.tmp ${gene1}-${gene2}-exonMappings-$paper.txt
	rm -rf oeg1.tmp oeg2.tmp oepairs.tmp int1.tmp int2.tmp int1.sorted.tmp int2.sorted.tmp 
done




exit


# extract from ExTraMapper results
#cat output/human-mouse/$gene1-$gene2/exonLevelMappings-1.0.txt | awk '{print $5"-"$11}' | sort > mypairs.tmp

