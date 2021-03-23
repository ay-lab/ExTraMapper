## Steps to generate the input files
The users should run the _run_human_mouse_preprocess_data_steps.sh_ to generate the inputfiles. All the input files will be generated under _preprocess/data_ folder. All the required executables and scripts are provided in the same folder. The _run_human_mouse_preprocess_data_steps.sh_ has six individual steps and should be run in the following manner 

### Run the following steps 

 - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) Set the following path <br>
   
    ```bash
    export EXTRAMAPPER_DIR=/path/to/Human-Mouse-Preprocess-Data/folder
    cd $EXTRAMAPPER_DIR
    ```
   <br>
 - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) Run the following step to fetch organism specific chromosomal fasta, gtf and liftOver files. <br>
   
    ```batch
    $ ./run_human_mouse_preprocess_data_steps.sh 0
    ```
    <br>
 - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The next step will create the genomedata object files. This step requires genomedata package which can be installed by running the following commnand. <br>
    
    ```bash
    $ pip install genomedata --user
    $ ./run_human_mouse_preprocess_data_steps.sh 1
    ```
    <br>
 - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The step below will create pickle files and gene summaries. The users are requested to install the latest pickle library. <br>
    
    ```bash
    $ ./run_human_mouse_preprocess_data_steps.sh 2
    ```
    <br>
 - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The following step will run the liftOver with multiple mappings. <br>
    
    ```bash
    $ ./run_human_mouse_preprocess_data_steps.sh 3
    ```
    <br>
 - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The next three steps will generate the input files. <br>
    
    ```bash
    $ ./run_human_mouse_preprocess_data_steps.sh 4
    $ ./run_human_mouse_preprocess_data_steps.sh 5
    $ ./run_human_mouse_preprocess_data_steps.sh 6
    ```
    <br>

#### Once finished the run_human_mouse_preprocess_data_steps.sh script shoudld produce the _data_ folder with the following subfolders.<br>
(Note: The _phastConsScores_ directory under _data_ is not essential so it will be empty)  

```bash 
.
|-- bin
|   |-- liftOver
|   `-- liftover-withMultiples
|-- data
|   |-- human-mouse
|   |   |-- GTFsummaries
|   |   |-- ensemblDownloads
|   |   |-- genePairsSummary-one2one.txt
|   |   |-- genomedataArchives
|   |   |-- liftoverRelatedFiles
|   |   |-- perExonLiftoverCoords
|   |   |-- perGenePairPickledInfo
|   |   `-- phastConsScores
|   |-- liftover_chains
|   |   |-- hg38
|   |   `-- mm10
|   `-- reference_genomes
|       |-- hg38
|       `-- mm10
`-- scripts
    `-- splitExonsIntoIndividualFiles.py
```
<br>

##### The whole process should take some time to finish!
##### [(Check also the Human-Moneky data processing steps)](https://github.com/ay-lab/ExTraMapper/tree/master/Human-Monkey-Processed-Data)
