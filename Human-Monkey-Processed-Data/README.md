## Steps to generate the input files
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
