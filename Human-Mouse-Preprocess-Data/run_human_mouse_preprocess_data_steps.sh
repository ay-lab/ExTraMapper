#!/bin/bash -ex
set -o pipefail
set -o errexit 

# read all the configuration parameters first
source config.conf
export EXTRAMAPPER_DIR=$EXTRAMAPPERDIR

step=$1

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

# Download necessary files from various resources
if [[ $step -eq 0 ]]; then
# download "per organism" specific files and keep the original organism names for future reuse
	## download the two reference genomes from UCSC and get rid of unknown, random and alt contigs
	for ref in $ref1 $ref2; do
		mkdir -p $referenceGenomesDir/$ref
		wget --timestamping ftp://hgdownload.cse.ucsc.edu/goldenPath/$ref/chromosomes/* --directory-prefix=$referenceGenomesDir/$ref
		rm -rf $referenceGenomesDir/$ref/chrUn* $referenceGenomesDir/$ref/*random* $referenceGenomesDir/$ref/*_alt*
		echo "Reference genome $ref downloaded to: $referenceGenomesDir/$ref"
	done

	# download liftover chains for each genome
	mkdir -p $chainsDir/$ref1/liftOver/ $chainsDir/$ref2/liftOver/

	## Assumes python version is 3 or greater ##
	ref1Cap=`echo $ref1 | python -c "s=input(); print (s[0].upper()+s[1:])"`
	ref2Cap=`echo $ref2 | python -c "s=input(); print (s[0].upper()+s[1:])"`
	wget http://hgdownload.cse.ucsc.edu/goldenPath/$ref1/liftOver/${ref1}To${ref2Cap}.over.chain.gz --directory-prefix=$chainsDir/$ref1/liftOver
	wget http://hgdownload.cse.ucsc.edu/goldenPath/$ref2/liftOver/${ref2}To${ref1Cap}.over.chain.gz --directory-prefix=$chainsDir/$ref2/liftOver
	chain1to2=`ls $chainsDir/${ref1}/liftOver/*.over.chain.gz`
	chain2to1=`ls $chainsDir/${ref2}/liftOver/*.over.chain.gz`
	echo "Liftover chain downloaded to: $chain1to2"
	echo "Liftover chain downloaded to: $chain2to1"

# download "per organism pair" files and name them org1 org2 to be generic
	# download the relevant Ensembl files for the pair of organisms of interest
	mkdir -p $ensemblDir
	wget ftp://ftp.ensembl.org/pub/release-$releaseNo/gtf/$org1EnsemblName/*.$releaseNo.gtf.gz -O $ensemblDir/org1.gtf.gz
	wget ftp://ftp.ensembl.org/pub/release-$releaseNo/gtf/$org2EnsemblName/*.$releaseNo.gtf.gz -O $ensemblDir/org2.gtf.gz
	echo "GTF files downloaded to: $ensemblDir"

	## You may get an error here, check org2EnsemblMartName or org2EnsemblMartNameShort applies ##
	wget ftp://ftp.ensembl.org/pub/release-$releaseNo/mysql/ensembl_mart_$releaseNo/${org1EnsemblMartName}_gene_ensembl__homolog_${org2EnsemblMartName}__dm.txt.gz -O $ensemblDir/org1_homolog_org2.txt.gz
	wget ftp://ftp.ensembl.org/pub/release-$releaseNo/mysql/ensembl_mart_$releaseNo/${org2EnsemblMartName}_gene_ensembl__homolog_${org1EnsemblMartName}__dm.txt.gz -O $ensemblDir/org2_homolog_org1.txt.gz
	echo "Ensembl homolog files downloaded to: $ensemblDir"


# Initialize the genomedata archives with the genomes of org1 and org2 
# Then add the phastConsScores to these archives
# This IS NOT dependent on Ensembl release number as long as the reference genome version DOES NOT change!
# Make sure genomedata is installed first 
elif [[ $step -eq 1 ]]; then

	##########################################################################################
	##### if genomedata package is not installed please do so using below command ##############
	##### for more info on genomedata: https://www.pmgenomics.ca/hoffmanlab/proj/genomedata/ ###
	# pip install genomedata --user
	############################################################################################

	mkdir -p $genomedataDir
	cd $genomedataDir 

	i=1
	for ref in $ref1 $ref2; do
		refdir=$referenceGenomesDir/$ref
		# Make sure there is no existing archive.
		rm -rf $ref $ref.fa
		#create a genome file with canonical names
		genome=$ref.fa
		for f in `ls $refdir/*.fa.gz`; do
			zcat $f >> $genome
		done
		# Make a genomedata archive containing just the sequence.
		genomedata-load-seq org$i $ref.fa
		genomedata-close-data org$i
		rm -rf $ref.fa
		i=$(($i+1))
	done

	##########################################################################################
	### This assumes the mafFiles are already generated locally #########################################
	#localPhastConsDir=/home/abhijit/overflow/proj_overflow/ExtraMapper/Python3_Version/automated_data_download/localphastconscores
	#for ref in $ref1 $ref2; do
	#	mkdir -p $phastConsDir/$ref-mafFiles/SCORES
	#	cp -r $localPhastConsDir/$ref-mafFiles/SCORES/* $phastConsDir/$ref-mafFiles/SCORES
	#done
	######################################################################################################
	
	#i=1
	#for ref in $ref1 $ref2; do
	#	wigsDir=$phastConsDir/$ref-mafFiles/SCORES
	#	genomedata-open-data org$i --tracknames phastCons
	#	for f in `ls $wigsDir/*.wig.gz`; do
	#		#echo $f
	#		zcat $f | genomedata-load-data org$i phastCons
	#	done
	#	genomedata-close-data org$i
	#	i=$(($i+1))
	#done
	cd -

# OPTIONAL step 2.1 (21 as int)
# Extract sequences of each exon as a fasta entry using genomedata archive
# Use whether the splice site signature is there to make sure coordinates are
# retrieved correctly for the fasta files
#elif [[ $step -eq 21 ]]; then
#	
#	##########################################################################################
#	##### if weblogo package is not installed please do so using below commands ##############
#	##### for more info on weblogo: http://weblogo.berkeley.edu/ #############################
#	#wget http://weblogo.berkeley.edu/release/weblogo.2.8.2.tar.gz
#	#tar -xzvf weblogo.2.8.2.tar.gz; rm -rf weblogo.2.8.2.tar.gz
#	##########################################################################################
#
#	seqlogobin=${EXTRAMAPPER_DIR}/preprocess/weblogo/seqlogo
#
#	outdir=${EXTRAMAPPER_DIR}/preprocess/output/seqLogos
#	python ${EXTRAMAPPER_DIR}/scripts/extractFastaSeqsForExons.py $dataDirPerPair $outdir
#
#	$seqlogobin -F PDF -f <(cat $outdir/org1_exonIntron_10-5p_10-3p.fasta | awk '{l1=$0; getline; printf("%s\t%s\n",l1,$0)}' \
#		| grep plusStrand |  awk '{print $1; print $2}')  -a -b -c -n -Y -E > $outdir/org1_exonIntron_10-5p_10-3p-onlyPlusStrand.pdf
#
#	$seqlogobin -F PDF -f <(cat $outdir/org1_exonIntron_10-5p_10-3p.fasta)  -a -b -c -n -Y -E > $outdir/org1_exonIntron_10-5p_10-3p-all.pdf
#	$seqlogobin -F PDF -f <(cat $outdir/org1_intronExon_10-5p_10-3p.fasta)  -a -b -c -n -Y -E > $outdir/org1_intronExon_10-5p_10-3p-all.pdf
#
#	$seqlogobin -F PDF -f <(cat $outdir/org2_exonIntron_10-5p_10-3p.fasta)  -a -b -c -n -Y -E > $outdir/org2_exonIntron_10-5p_10-3p-all.pdf
#	$seqlogobin -F PDF -f <(cat $outdir/org2_intronExon_10-5p_10-3p.fasta)  -a -b -c -n -Y -E > $outdir/org2_intronExon_10-5p_10-3p-all.pdf
#
# extract conservation scores and acceptor-donor sites for each exon
elif [[ $step -eq 2 ]]; then

        ## Extract the gz files ##
	gunzip -k $ensemblDir/org1.gtf.gz
	gunzip -k $ensemblDir/org2.gtf.gz
	gunzip -k $ensemblDir/org1_homolog_org2.txt.gz
	gunzip -k $ensemblDir/org2_homolog_org1.txt.gz	

        python ${EXTRAMAPPER_DIR}/scripts/parseAndPicklePerPair.py $dataDirPerPair $GTFsummaryDir $perGenePairPickleDir
	mv $perGenePairPickleDir/genePairsSummary-one2one.txt $dataDirPerPair/genePairsSummary-one2one.txt


# liftOver the exon lists but this time allow multiple mappings and also compute intersections with the other set of exons
elif [[ $step -eq 3 ]]; then
	# liftover either all exons or just the ones we cared about
	indir=$GTFsummaryDir # all exons
	# OR
	indir=$GTFsummaryDir/onlyOrthologAndCodingGenes # ~15.8k ortho and protein coding gene pairs

	mkdir -p $liftOverFilesDir
	# extract the bed files for exons (all, nonCoding, partCoding, and fullCoding)
	cat $indir/org1-allExons-GTFparsed.txt | awk 'NR>1{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | sort -k1,1 -k2,2n > $liftOverFilesDir/org1_allExonsList.bed
	cat $indir/org2-allExons-GTFparsed.txt | awk 'NR>1{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | sort -k1,1 -k2,2n > $liftOverFilesDir/org2_allExonsList.bed
	cat $indir/org1-allExons-GTFparsed.txt |\
		awk '$6=="partCoding" {print $1"\t"$7"\t"$8"\t"$4"\t"$5}' | sort -k1,1 -k2,2n >$liftOverFilesDir/org1_partCodingExonsList.bed
	cat $indir/org2-allExons-GTFparsed.txt |\
		awk '$6=="partCoding" {print $1"\t"$7"\t"$8"\t"$4"\t"$5}' | sort -k1,1 -k2,2n >$liftOverFilesDir/org2_partCodingExonsList.bed

	cat $indir/org1-allExons-GTFparsed.txt | awk '$6=="fullCoding" {print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > $liftOverFilesDir/org1_f.temp
	cat $indir/org2-allExons-GTFparsed.txt | awk '$6=="fullCoding" {print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > $liftOverFilesDir/org2_f.temp

	cat $liftOverFilesDir/org1_partCodingExonsList.bed $liftOverFilesDir/org1_f.temp | sort -k1,1 -k2,2n > $liftOverFilesDir/org1_allCodingExonsList.bed
	cat $liftOverFilesDir/org2_partCodingExonsList.bed $liftOverFilesDir/org2_f.temp | sort -k1,1 -k2,2n > $liftOverFilesDir/org2_allCodingExonsList.bed

	cat $liftOverFilesDir/org1_allCodingExonsList.bed | awk '{print $5,$0}' | sort -k1,1 > $liftOverFilesDir/org1_allCodingExonsList.sorted.temp
	cat $liftOverFilesDir/org2_allCodingExonsList.bed | awk '{print $5,$0}' | sort -k1,1 > $liftOverFilesDir/org2_allCodingExonsList.sorted.temp
	cat $liftOverFilesDir/org1_allExonsList.bed | awk '{print $5,$0}' | sort -k1,1 > $liftOverFilesDir/org1_allExonsList.sorted.temp
	cat $liftOverFilesDir/org2_allExonsList.bed | awk '{print $5,$0}' | sort -k1,1 > $liftOverFilesDir/org2_allExonsList.sorted.temp
	cat $liftOverFilesDir/org1_partCodingExonsList.bed | awk '{print $5,$0}' | sort -k1,1 > $liftOverFilesDir/org1_partCodingExonsList.sorted.temp
	cat $liftOverFilesDir/org2_partCodingExonsList.bed | awk '{print $5,$0}' | sort -k1,1 > $liftOverFilesDir/org2_partCodingExonsList.sorted.temp
	
	chain1to2=`ls $chainsDir/${ref1}/liftOver/*.over.chain.gz`
        chain2to1=`ls $chainsDir/${ref2}/liftOver/*.over.chain.gz`

	for flank in 0; do 
		for minMatch in 1.0 0.95 0.9; do
			${EXTRAMAPPER_DIR}/preprocess/bin/liftover-withMultiples $flank $minMatch $chain1to2 $chain2to1
		done
	done
	rm -rf $liftOverFilesDir/org2_allExonsList.sorted.temp $liftOverFilesDir/org1_allExonsList.sorted.temp 
	rm -rf $liftOverFilesDir/org2_partCodingExonsList.sorted.temp $liftOverFilesDir/org1_partCodingExonsList.sorted.temp
	rm -rf $liftOverFilesDir/org2_allCodingExonsList.sorted.temp $liftOverFilesDir/org1_allCodingExonsList.sorted.temp

# put together, sort, uniq and then split into one file per exon all the liftover files created so far
elif [[ $step -eq 4 ]]; then
	indir=$liftOverFilesDir
	outdir=$perExonLiftoverDir
	flank=0
	
	rm -rf oneHugeFile-2to1-partCoding.txt oneHugeFile-1to2-partCoding.txt
	for minMatch in 1.0 0.95 0.9; do
		suffix=flank$flank-minMatch$minMatch-multiples-partCoding
		zcat $indir/org1_VS_org2_to_org1_intersectingExonsList-$suffix.bed.gz  |\
			awk '$6!="."{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$11"\t"$9"\t"$12"\t"s}' s=$suffix >> oneHugeFile-2to1-partCoding.txt
		zcat $indir/org2_VS_org1_to_org2_intersectingExonsList-$suffix.bed.gz  |\
			awk '$6!="."{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$11"\t"$9"\t"$12"\t"s}' s=$suffix >> oneHugeFile-1to2-partCoding.txt
	done

	rm -rf oneHugeFile-2to1-others.txt oneHugeFile-1to2-others.txt
	for minMatch in 1.0 0.95 0.9; do
		suffix=flank$flank-minMatch$minMatch-multiples
		zcat $indir/org1_VS_org2_to_org1_intersectingExonsList-$suffix.bed.gz  |\
			awk '$6!="."{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$11"\t"$9"\t"$12"\t"s}' s=$suffix >> oneHugeFile-2to1-others.txt
		zcat $indir/org2_VS_org1_to_org2_intersectingExonsList-$suffix.bed.gz  |\
			awk '$6!="."{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$11"\t"$9"\t"$12"\t"s}' s=$suffix >> oneHugeFile-1to2-others.txt
	done

	#cat oneHugeFile-1to2-partCoding.txt | sort -u -k1,10 | sort -k10,10 > oneHugeFile-1to2-partCoding.txt.sorted
	#cat oneHugeFile-2to1-partCoding.txt | sort -u -k1,10 | sort -k10,10 > oneHugeFile-2to1-partCoding.txt.sorted
	#rm -rf oneHugeFile-1to2-partCoding.txt oneHugeFile-2to1-partCoding.txt

	#cat oneHugeFile-1to2-others.txt | sort -u -k1,10 | sort -k10,10 > oneHugeFile-1to2-others.txt.sorted
	#cat oneHugeFile-2to1-others.txt | sort -u -k1,10 | sort -k10,10 > oneHugeFile-2to1-others.txt.sorted
	#rm -rf oneHugeFile-1to2-others.txt oneHugeFile-2to1-others.txt

	cat oneHugeFile-1to2-partCoding.txt oneHugeFile-1to2-others.txt | sort -k10,10 >oneHugeFile-1to2.txt.sorted	
	cat oneHugeFile-2to1-partCoding.txt oneHugeFile-2to1-others.txt | sort -k10,10 >oneHugeFile-2to1.txt.sorted	

	mkdir -p $outdir/org1 $outdir/org2
	whichCol=10 # 10th column will be 9th in 0-based python indices
	fileSuffix="_mapped.txt"
	python ${EXTRAMAPPER_DIR}/preprocess/scripts/splitExonsIntoIndividualFiles.py oneHugeFile-1to2.txt.sorted $outdir/org1 $whichCol $fileSuffix
	python ${EXTRAMAPPER_DIR}/preprocess/scripts/splitExonsIntoIndividualFiles.py oneHugeFile-2to1.txt.sorted $outdir/org2 $whichCol $fileSuffix
	rm -rf oneHugeFile*.txt

	## NOTE THAT BELOW SPLIT BY AWK WON'T WORK FOR SUCH MANY FILES, AWK KEEPS HANDLES OPEN!
	#mkdir -p exonLiftoverCoords/human exonLiftoverCoords/mouse
	#cat oneHugeFile-2to1.txt | awk '{print > d"/"$10".txt"}' d="exonLiftoverCoords/mouse"
	#cat oneHugeFile-1to2.txt | awk '{print > d"/"$10".txt"}' d="exonLiftoverCoords/human"

# put together, sort, uniq and then split into one file per exon all the liftover files for UNMAPPED EXONS so far
elif [[ $step -eq 5 ]]; then
	indir=$liftOverFilesDir
	outdir=$perExonLiftoverDir
	flank=0

	rm -rf oneHugeFile-2to1-partCoding.txt oneHugeFile-1to2-partCoding.txt
	for minMatch in 1.0 0.95 0.9; do
		suffix=flank$flank-minMatch$minMatch-multiples-partCoding
		zcat $indir/org1_to_org2_liftOver_unmappedExonsList-$suffix.bed.gz |\
			awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"s}' s=$suffix >> oneHugeFile-1to2-partCoding.txt
		zcat $indir/org2_to_org1_liftOver_unmappedExonsList-$suffix.bed.gz |\
			awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"s}' s=$suffix >> oneHugeFile-2to1-partCoding.txt
	done

	rm -rf oneHugeFile-2to1-others.txt oneHugeFile-1to2-others.txt
	for minMatch in 1.0 0.95 0.9; do
		suffix=flank$flank-minMatch$minMatch-multiples
		zcat $indir/org1_to_org2_liftOver_unmappedExonsList-$suffix.bed.gz |\
			awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"s}' s=$suffix >> oneHugeFile-1to2-others.txt
		zcat $indir/org2_to_org1_liftOver_unmappedExonsList-$suffix.bed.gz |\
			awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"s}' s=$suffix >> oneHugeFile-2to1-others.txt
	done

	#sort on the same field that goes to whichCol below
	cat oneHugeFile-1to2-partCoding.txt oneHugeFile-1to2-others.txt | sort -k5,5 >oneHugeFile-1to2.txt.sorted	
	cat oneHugeFile-2to1-partCoding.txt oneHugeFile-2to1-others.txt | sort -k5,5 >oneHugeFile-2to1.txt.sorted	

	mkdir -p $outdir/org1 $outdir/org2
	whichCol=5 # 10th column will be 9th in 0-based python indices
	fileSuffix="_unmapped.txt"
	python ${EXTRAMAPPER_DIR}/preprocess/scripts/splitExonsIntoIndividualFiles.py oneHugeFile-1to2.txt.sorted $outdir/org1 $whichCol $fileSuffix
	python ${EXTRAMAPPER_DIR}/preprocess/scripts/splitExonsIntoIndividualFiles.py oneHugeFile-2to1.txt.sorted $outdir/org2 $whichCol $fileSuffix

# put together, sort, uniq and then split into one file per exon all the liftover files for MAPPED EXONS that DO NOT INTERSECT WITH AN EXON so far
elif [[ $step -eq 6 ]]; then
	indir=$liftOverFilesDir
	outdir=$perExonLiftoverDir
	flank=0

	## !!!!!!!!!!! don't do the part coding versions as they will lose the splice site info once liftedover
	#rm -rf oneHugeFile-2to1-partCoding.txt oneHugeFile-1to2-partCoding.txt
	#for minMatch in 1.0 0.95 0.9; do
	#	suffix=flank$flank-minMatch$minMatch-multiples-partCoding
	#	zcat $indir/org1_VS_org2_to_org1_nonintersectingExonsList-$suffix.bed.gz |\
	#		awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"s}' s=$suffix >> oneHugeFile-2to1-partCoding.txt
	#	zcat $indir/org2_VS_org1_to_org2_nonintersectingExonsList-$suffix.bed.gz |\
	#		awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"s}' s=$suffix >> oneHugeFile-1to2-partCoding.txt
	#done

	rm -rf oneHugeFile-2to1-others.txt oneHugeFile-1to2-others.txt
	for minMatch in 1.0 0.95 0.9; do
		suffix=flank$flank-minMatch$minMatch-multiples
		zcat $indir/org1_VS_org2_to_org1_nonintersectingExonsList-$suffix.bed.gz |\
			awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"s}' s=$suffix >> oneHugeFile-2to1-others.txt
		zcat $indir/org2_VS_org1_to_org2_nonintersectingExonsList-$suffix.bed.gz |\
			awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"s}' s=$suffix >> oneHugeFile-1to2-others.txt
	done

	##sort on the same field that goes to whichCol below
	#cat oneHugeFile-1to2-partCoding.txt oneHugeFile-1to2-others.txt | sort -k5,5 >oneHugeFile-1to2.txt.sorted	
	#cat oneHugeFile-2to1-partCoding.txt oneHugeFile-2to1-others.txt | sort -k5,5 >oneHugeFile-2to1.txt.sorted	
	cat oneHugeFile-1to2-others.txt | sort -k5,5 >oneHugeFile-1to2.txt.sorted	
	cat oneHugeFile-2to1-others.txt | sort -k5,5 >oneHugeFile-2to1.txt.sorted	

	mkdir -p $outdir/org1 $outdir/org2
	whichCol=5 # 10th column will be 9th in 0-based python indices
	fileSuffix="_nonintersecting.txt"
	python ${EXTRAMAPPER_DIR}/preprocess/scripts/splitExonsIntoIndividualFiles.py oneHugeFile-1to2.txt.sorted $outdir/org1 $whichCol $fileSuffix
	python ${EXTRAMAPPER_DIR}/preprocess/scripts/splitExonsIntoIndividualFiles.py oneHugeFile-2to1.txt.sorted $outdir/org2 $whichCol $fileSuffix

	rm -rf oneHugeFile* dummy.txt
fi

exit
