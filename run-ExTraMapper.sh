#!/bin/bash -ex
set -o pipefail
set -o errexit 

# read all the configuration parameters first
source config.conf
export EXTRAMAPPER_DIR=$EXTRAMAPPERDIR

#define the variables
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

# Some example input gene pairs #
genePairId=ENSG00000078900-ENSMUSG00000029026 # TP73
genePairId=ENSG00000049768-ENSMUSG00000039521 # FOXP3
genePairId=ENSG00000157764-ENSMUSG00000002413 # BRAF
genePairId=ENSG00000143315-ENSMUSG00000050229 # 1 exon gene pair
genePairId=ENSG00000140416-ENSMUSG00000032366 # TPM1 vs Tpm1 
genePairId=ENSG00000141510-ENSMUSG00000059552 # TP53
#genePairId=$1

#define the liftover minimum match threshold
MAPPED_EXON_THRES=1.0

#run ExTramapper (python 2.7 tested)
python ${EXTRAMAPPER_DIR}/scripts/ExTraMapper.py ${MAPPED_EXON_THRES} $org1 $org2 $genePairId > allOutputs.txt
