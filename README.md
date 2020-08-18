# ExTraMapper
ExTraMapper is a tool to find Exon and Transcript-level Mappings of a given pair of orthologous genes between two organisms using sequence conservation. The figure below shows the overall schematic description of ExTraMapper mapping the homologous transcript and exon-pairs between human and mouse genome. 

![ExTraMapper_Figure](https://user-images.githubusercontent.com/18036388/90572310-8b693e00-e168-11ea-9fbc-8188c2834de9.jpg)

# Run ExtraMapper
### Read all the configuration parameters first
source config.conf

export EXTRAMAPPER_DIR=$EXTRAMAPPERDIR

### Set the liftover minimum match threshold
MAPPED_EXON_THRES=1.0

### An example input gene pair (BRAF)
genePairId="ENSG00000157764-ENSMUSG00000002413"

### Run ExTramapper (Requires python 2.7)
python $EXTRAMAPPER_DIR/scripts/ExTraMapper.py $MAPPED_EXON_THRES $org1 $org2 $genePairId > Outputs.txt

### Check the run-ExTraMapper.sh for more details
