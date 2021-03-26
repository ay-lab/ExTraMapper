## Steps to generate the input files (Human-Rhesus_macaque)
The users should run the _run_human_monkey_preprocess_data_steps.sh_ to generate the inputfiles. All the input files will be generated under _preprocess/data_ folder. All the required executables and scripts are provided in the same folder. The _run_human_monkey_preprocess_data_steps.sh_ has six individual steps and should be run in the following manner 

### Run the following steps 

 - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) Set the following path <br>
   
    ```bash
    export EXTRAMAPPER_DIR=/path/to/Human-Monkey-Processed-Data/folder
    cd $EXTRAMAPPER_DIR
    ```
   <br>
   
 - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) Run the following step to fetch organism specific chromosomal fasta, gtf and liftOver files. <br>
   
    ```batch
    $ ./run_human_monkey_preprocess_data_steps.sh 0
    ```
    <br>
    
 - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The above step may produce an error while downloading the monkey genome from UCSC. In that case, please do the following and the script will produce the required fasta files. <br>
    
    ```bash
    $ cd ./preprocess/data/reference_genomes/rheMac10/
    $ perl getFasta.pl
    $ cd -
    ```
    <br>
    
 - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The next step will create the genomedata object files. This step requires genomedata package which can be installed by running the following commnand. <br>
    
    ```bash
    $ pip install genomedata --user
    $ ./run_human_monkey_preprocess_data_steps.sh 1
    ```
    <br>
    
 - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The step below will create pickle files and gene summaries. The users are requested to install the latest pickle library. <br>
    
    ```bash
    $ ./run_human_monkey_preprocess_data_steps.sh 2
    ```
    <br>
    
 - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The following step will run the liftOver with multiple mappings. <br>
    
    ```bash
    $ ./run_human_monkey_preprocess_data_steps.sh 3
    ```
    <br>
    
 - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The next three steps will generate the input files. <br>
    
    ```bash
    $ ./run_human_mouse_preprocess_data_steps.sh 4
    $ ./run_human_mouse_preprocess_data_steps.sh 5
    $ ./run_human_mouse_preprocess_data_steps.sh 6
    ```
    <br>
    
#### Once finished the run_human_monkey_preprocess_data_steps.sh script shoudld produce the _data_ folder with the following subfolders.<br>

```bash 
.
|-- bin
|   |-- liftOver
|   `-- liftover-withMultiples
|-- data
|   |-- human-rhesus
|   |   |-- GTFsummaries
|   |   |   |-- onlyOrthologAndCodingGenes
|   |   |   |   |-- org1-allExons-GTFparsed.txt
|   |   |   |   |-- org1-allGenes-GTFparsed.txt
|   |   |   |   |-- org1-allTranscripts-GTFparsed.txt
|   |   |   |   |-- org2-allExons-GTFparsed.txt
|   |   |   |   |-- org2-allGenes-GTFparsed.txt
|   |   |   |   `-- org2-allTranscripts-GTFparsed.txt
|   |   |   |-- org1-allExons-GTFparsed.txt
|   |   |   |-- org1-allGenes-GTFparsed.txt
|   |   |   |-- org1-allTranscripts-GTFparsed.txt
|   |   |   |-- org2-allExons-GTFparsed.txt
|   |   |   |-- org2-allGenes-GTFparsed.txt
|   |   |   `-- org2-allTranscripts-GTFparsed.txt
|   |   |-- ensemblDownloads
|   |   |   |-- org1.gtf
|   |   |   |-- org1.gtf.gz
|   |   |   |-- org1_homolog_org2.txt
|   |   |   |-- org1_homolog_org2.txt.gz
|   |   |   |-- org2.gtf
|   |   |   |-- org2.gtf.gz
|   |   |   |-- org2_homolog_org1.txt
|   |   |   `-- org2_homolog_org1.txt.gz
|   |   |-- genePairsSummary-one2one.txt
|   |   |-- genomedataArchives
|   |   |   |-- org1 [25 entries exceeds filelimit, not opening dir]
|   |   |   `-- org2 [22 entries exceeds filelimit, not opening dir]
|   |   |-- liftoverRelatedFiles [56 entries exceeds filelimit, not opening dir]
|   |   |-- perExonLiftoverCoords
|   |   |   |-- org1 [619127 entries exceeds filelimit, not opening dir]
|   |   |   `-- org2 [260616 entries exceeds filelimit, not opening dir]
|   |   `-- perGenePairPickledInfo [16150 entries exceeds filelimit, not opening dir]
|   |-- liftover_chains
|   |   |-- hg38
|   |   |   `-- liftOver
|   |   |       `-- hg38ToRheMac10.over.chain.gz
|   |   `-- rheMac10
|   |       `-- liftOver
|   |           `-- rheMac10ToHg38.over.chain.gz
|   `-- reference_genomes
|       |-- hg38 [27 entries exceeds filelimit, not opening dir]
|       `-- rheMac10 [25 entries exceeds filelimit, not opening dir]
`-- scripts
    `-- splitExonsIntoIndividualFiles.py
```

##### The whole process should take some time to finish!
##### [(Check also the Human-Mouse data processing steps)](https://github.com/ay-lab/ExTraMapper/tree/master/Human-Mouse-Preprocess-Data)
