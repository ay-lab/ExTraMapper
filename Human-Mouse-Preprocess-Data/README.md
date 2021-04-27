## Steps to generate the input files (Human - Mouse)
The users should run the _extMpreprocess_ to generate the inputfiles. All the input files will be generated under _preprocess/data_ folder. All the required executables and scripts are provided here. The _extMpreprocess_ has 7 individual steps and should be run in the following manner 

### Run the following steps 

 - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) For help, type <br>
   
    ```bash
    ./extMpreprocess help
    
    This script will download and preprocess the dataset required for exon-pair and transcript pair finding by ExTraMapper.
    Type ./extMpreprocess <config.conf> <step> to execute the script.
    Type ./extMpreprocess example to print a example config.conf file.

    This script will run seven (7) sequential steps to create the inputs for ExTraMapper program.
    Users can provide step numbers (1-7) or all in the <step> arugemt of this script.
    Short description of the individual scripts:
    Step 1: Download per organism specific files e.g. reference genomes, gene annotation files.
    Step 2: Will create genomedata archives with the genomes of org1 and org2 (Make sure to install genomedata package).
    Step 3: Pickle files for each homologous gene pair will be created.
    Step 4: Perform coordinate liftOver of exons with multiple mappings (This step requires bedtools and liftOver executables).
    Step 5-7: postprocessing the liftOver files.
    ```
   <br>
   <br>
 - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The script requires genomedata package which can be installed by running the following commnand. <br>
    
    ```bash
    $ pip install genomedata --user
    $ ./run_human_mouse_preprocess_data_steps.sh 1
    ```
    <br>
 
    <br>

#### Once finished the run_human_mouse_preprocess_data_steps.sh script shoudld produce the _data_ folder with the following subfolders.<br>

```bash 
./preprocess
|-- bin
|   `-- liftOver
|-- data
    |-- human-mouse
    |   |-- GTFsummaries
    |   |   |-- onlyOrthologAndCodingGenes
    |   |   |   |-- org1-allExons-GTFparsed.txt
    |   |   |   |-- org1-allGenes-GTFparsed.txt
    |   |   |   |-- org1-allTranscripts-GTFparsed.txt
    |   |   |   |-- org2-allExons-GTFparsed.txt
    |   |   |   |-- org2-allGenes-GTFparsed.txt
    |   |   |   `-- org2-allTranscripts-GTFparsed.txt
    |   |   |-- org1-allExons-GTFparsed.txt
    |   |   |-- org1-allGenes-GTFparsed.txt
    |   |   |-- org1-allTranscripts-GTFparsed.txt
    |   |   |-- org2-allExons-GTFparsed.txt
    |   |   |-- org2-allGenes-GTFparsed.txt
    |   |   `-- org2-allTranscripts-GTFparsed.txt
    |   |-- ensemblDownloads
    |   |   |-- org1.gtf
    |   |   |-- org1.gtf.gz
    |   |   |-- org1_homolog_org2.txt
    |   |   |-- org1_homolog_org2.txt.gz
    |   |   |-- org2.gtf
    |   |   |-- org2.gtf.gz
    |   |   |-- org2_homolog_org1.txt
    |   |   `-- org2_homolog_org1.txt.gz
    |   |-- genePairsSummary-one2one.txt
    |   |-- genomedataArchives
    |   |   |-- org1 [25 entries exceeds filelimit, not opening dir]
    |   |   `-- org2 [22 entries exceeds filelimit, not opening dir]
    |   |-- liftoverRelatedFiles [56 entries exceeds filelimit, not opening dir]
    |   |-- perExonLiftoverCoords
    |   |   |-- org1 [654707 entries exceeds filelimit, not opening dir]
    |   |   `-- org2 [484860 entries exceeds filelimit, not opening dir]
    |   |-- perGenePairPickledInfo [15804 entries exceeds filelimit, not opening dir]
    |   
    |-- liftover_chains
    |   |-- hg38
    |   |   `-- liftOver
    |   |       `-- hg38ToMm10.over.chain.gz
    |   `-- mm10
    |       `-- liftOver
    |           `-- mm10ToHg38.over.chain.gz
    `-- reference_genomes
        |-- hg38 [27 entries exceeds filelimit, not opening dir]
        `-- mm10 [24 entries exceeds filelimit, not opening dir]

```
<br>

##### The whole process should take some time to finish!
##### [(Check also the Human-Moneky data processing steps)](https://github.com/ay-lab/ExTraMapper/tree/master/Human-Monkey-Processed-Data)
