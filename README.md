# ExTraMapper
ExTraMapper is a tool to find Exon and Transcript-level Mappings of a given pair of orthologous genes between two organisms using sequence conservation. The figure below shows the overall schematic description of ExTraMapper mapping the homologous transcript and exon-pairs between human and mouse genome. 


![ExTraMapper_Figure](https://user-images.githubusercontent.com/18036388/90572310-8b693e00-e168-11ea-9fbc-8188c2834de9.jpg)

# Steps to run ExtraMapper (For python version 3 or later usage)

### Step 1: Prepare the input files
ExTraMapper requires a set of preprocessed files to find the conservation scores. Examples to create these files are provided within the following folders
 
1. [__Human-Mouse-Preprocessed-Data__](https://github.com/ay-lab/ExTraMapper/tree/master/Human-Mouse-Processed-Data) 

    Quick look:
    
   - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) Set the following path <br><br>
   
    ```bash
    export EXTRAMAPPER_DIR=/path/to/Human-Mouse-Preprocess-Data/folder
    cd $EXTRAMAPPER_DIR
    ```
   
    - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) Run the following step to fetch organism specific chromosomal fasta, gtf and liftOver files. <br><br>
   
    ```batch
    $ ./run_human_mouse_preprocess_data_steps.sh 0
    ```
    
    - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The next step will create the genomedata object files. This step requires genomedata package which can be installed by running the following commnand. <br><br>
    
    ```bash
    $ pip install genomedata --user
    $ ./run_human_mouse_preprocess_data_steps.sh 1
    ```
    
    - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The step below will create pickle files and gene summaries. The users are requested to install the latest pickle library. <br><br>
    
    ```bash
    $ ./run_human_mouse_preprocess_data_steps.sh 2
    ```
    
    - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The following step will run the liftOver with multiple mappings. <br><br>
    
    ```bash
    $ ./run_human_mouse_preprocess_data_steps.sh 3
    ```
    
    - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The next three steps will generate the input files. <br><br>
    
    ```bash
    $ ./run_human_mouse_preprocess_data_steps.sh 4
    $ ./run_human_mouse_preprocess_data_steps.sh 5
    $ ./run_human_mouse_preprocess_data_steps.sh 6
    ```
    
    and 
    
2. [__Human-Mokey-Preprocessed-Data__](https://github.com/ay-lab/ExTraMapper/tree/master/Human-Monkey-Processed-Data) 

    Quick look:
    
   - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) Set the following path <br><br>
   
    ```bash
    export EXTRAMAPPER_DIR=/path/to/Human-Monkey-Preprocess-Data/folder
    cd $EXTRAMAPPER_DIR
    ```
   
    - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) Run the following step to fetch organism specific chromosomal fasta, gtf and liftOver files. <br><br>
   
    ```batch
    $ ./run_human_monkey_preprocess_data_steps.sh 0
    ```
    
    - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The above step may produce an error while downloading the monkey genome from UCSC. In that case, please do the following and the script will produce the required fasta files. <br><br>
    
    ```bash
    $ cd ./preprocess/data/reference_genomes/rheMac10/
    $ perl getFasta.pl
    $ cd -
    ```
    
    - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The next step will create the genomedata object files. This step requires genomedata package which can be installed by running the following commnand. <br><br>
    
    ```bash
    $ pip install genomedata --user
    $ ./run_human_monkey_preprocess_data_steps.sh 1
    ```
    
    - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The step below will create pickle files and gene summaries. The users are requested to install the latest pickle library. <br><br>
    
    ```bash
    $ ./run_human_monkey_preprocess_data_steps.sh 2
    ```
    
    - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The following step will run the liftOver with multiple mappings. <br><br>
    
    ```bash
    $ ./run_human_monkey_preprocess_data_steps.sh 3
    ```
    
    - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) The next three steps will generate the input files. <br><br>
    
    ```bash
    $ ./run_human_monkey_preprocess_data_steps.sh 4
    $ ./run_human_monkey_preprocess_data_steps.sh 5
    $ ./run_human_monkey_preprocess_data_steps.sh 6
    ```
   
Users should look into these folders and follow the instructions to create the required input files before going to the next step.   


### Step 2: Set the following path
```bash export EXTRAMAPPER_DIR=/path/to/this/folder```

### Step 3: Run ExTraMapper individually
```bash
$ python ExTraMapper.py -h
usage: ExTraMapper.py [-h] -m MAPPING -o1 ORG1 -o2 ORG2 -p ORTHOLOG

Check the help flag

optional arguments:
  -h, --help   show this help message and exit
  -m MAPPING   ExTraMapper Exon threshold value [e.g. 1]
  -o1 ORG1     First organism name [e.g. human]
  -o2 ORG2     Second organism name [e.g. mouse]
  -p ORTHOLOG  Orthologous gene pair [e.g. ENSG00000141510-ENSMUSG00000059552 OR all]
```

# OR

### Step 3: Run ExTraMapper for all the gene pairs
```bash
$ python ExTraMapper.py -h
usage: ExTraMapper.py [-h] -m MAPPING -o1 ORG1 -o2 ORG2 -p all
```

### Refer the work
[_ExTraMapper: Exon- and Transcript-level mappings for orthologous gene pairs._](https://www.biorxiv.org/content/10.1101/277723v1)

The data shown in the above paper was performed using Human & Mouse ENSMBL release 81 with python 2.7 code. 
The current update is with ENSMBL release 102 and python 3 or later version. To see the older code and data please
change the __Branch__ to __ExTraMapper-python2v__ from __master__
