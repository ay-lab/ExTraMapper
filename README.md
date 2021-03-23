# ExTraMapper
ExTraMapper is a tool to find Exon and Transcript-level Mappings of a given pair of orthologous genes between two organisms using sequence conservation. The figure below shows the overall schematic description of ExTraMapper mapping the homologous transcript and exon-pairs between human and mouse genome. 


![ExTraMapper_Figure](https://user-images.githubusercontent.com/18036388/90572310-8b693e00-e168-11ea-9fbc-8188c2834de9.jpg)

# Steps to run ExtraMapper (For python version 3 or later usage)

### Step 1: Prepare the input files
ExTraMapper requires a set of preprocessed files to find the conservation scores. Examples to create these files are provided within the following folders
 
1. [__Human-Mouse-Preprocessed-Data__](https://github.com/ay-lab/ExTraMapper/tree/master/Human-Mouse-Preprocess-Data) 

    Quick look:
    
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
    
    and 
    
2. [__Human-Mokey-Preprocessed-Data__](https://github.com/ay-lab/ExTraMapper/tree/master/Human-Monkey-Processed-Data) 

    Quick look:
    
   - ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) Set the following path <br>
   
    ```bash
    export EXTRAMAPPER_DIR=/path/to/Human-Monkey-Preprocess-Data/folder
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
    $ ./run_human_monkey_preprocess_data_steps.sh 4
    $ ./run_human_monkey_preprocess_data_steps.sh 5
    $ ./run_human_monkey_preprocess_data_steps.sh 6
    ```
   <br>
   
Users should look into these folders and follow the instructions to create the required input files before going to the next step.   

<br>

### Step 2: Set the following path
```bash export EXTRAMAPPER_DIR=/path/to/this/folder```

<br>

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

#### Example run of ExTraMapper.py using orthologous gene pair ENSG00000141510-ENSMUSG00000059552 
```bash
$ python ExTraMapper.py -m 1 -o1 human -o2 mouse -p ENSG00000141510-ENSMUSG00000059552

Finding exon mappings for gene pair number 0    ENSG00000141510-ENSMUSG00000059552
*****************************************************************
Gene pair ID: ENSG00000141510-ENSMUSG00000059552

Information about each gene. Last two numbers are no of transcripts and exons
ENSG00000141510 chr17   7661779 7687538 -       ENSG00000141510 TP53    protein_coding  27      49      gene
ENSMUSG00000059552      chr11   69580359        69591873        +       ENSMUSG00000059552      Trp53   protein_coding  6       24      gene

Number of exons before and after duplicate removal according to coordinates
Org1    49      40
Org2    24      20

*****************************************************************

*****************************************************************
GCGCTGGGGACCTGTCCCTAGGGGGCAGATGAGACACTGATGGGCGTACTTAGAGATTTGCCATGAAGTGGGTTTGAAGAATGGAGCTGTGTGTGAAAT
Exon file type summaries for the first gene from: ENSG00000141510-ENSMUSG00000059552
        0 exons with: No file exists
        22 exons with: Only Mapped
        0 exons with: Only nonintersecting
        11 exons with: Only unmapped
        15 exons with: Mapped and unmapped
        0 exons with: Mapped and nonintersecting
        1 exons with: Nonintersecting and unmapped
        0 exons with: All three files
Exon file type summaries for the second gene from: ENSG00000141510-ENSMUSG00000059552
        0 exons with: No file exists
        14 exons with: Only Mapped
        0 exons with: Only nonintersecting
        3 exons with: Only unmapped
        7 exons with: Mapped and unmapped
        0 exons with: Mapped and nonintersecting
        0 exons with: Nonintersecting and unmapped
        0 exons with: All three files
Writing exon-level similarity scores into file:
 /path/output/human-mouse/ENSG00000141510-ENSMUSG00000059552/exonLevelSimilarities-1.0.txt

Writing exon classes into file:
 /path/output/human-mouse/ENSG00000141510-ENSMUSG00000059552/exonClasses-1.0.txt
        For org1: Mapped exons= 17, Unmapped exons= 21, Nonintersecting exons= 1, OTHER= 10
        For org2: Mapped exons= 13, Unmapped exons= 7, Nonintersecting exons= 0, OTHER= 4
*****************************************************************

*****************************************************************
Writing exon-level mappings into file:
 /path/output/human-mouse/ENSG00000141510-ENSMUSG00000059552/exonLevelMappings-1.0.txt
Writing trascript-level similarity scores into file:
 /path/output/human-mouse/ENSG00000141510-ENSMUSG00000059552/transcriptLevelSimilarities-1.0.txt
Writing transcript-level mappings into file:
 /path/output/human-mouse/ENSG00000141510-ENSMUSG00000059552/transcriptLevelMappings-1.0.txt

Condition counter from the greedy transcript mapping stage:
        5 pairs with Condition1: Unique winner pair
        0 pairs with Condition2: Tie in one score, not in the other
        0 pairs with Condition3: Tie in both scores but coding exon length diff breaks the tie
        0 pairs with Condition4: Tie in both scores and coding exon length diff but overall exon length breaks the tie
        1 pairs with Condition5: Tie in all the above but coding length (bp) diff breaks the tie
        0 pairs with Condition6: Tie in all the above, just give up and report all

Writing UCSC browser bed output for org1 into file:
 /path/output/human-mouse/ENSG00000141510-ENSMUSG00000059552/org1-ucsc-1.0.bed
Writing UCSC browser bed output for org2 into file:
 /path/output/human-mouse/ENSG00000141510-ENSMUSG00000059552/org2-ucsc-1.0.bed

........
ExTraMapper ran successfully for 1 gene pairs between: human and mouse


*****************************************************************
$ tree ./output

./output
`-- human-mouse
    `-- ENSG00000141510-ENSMUSG00000059552
        |-- exonClasses-1.0.txt
        |-- exonLevelMappings-1.0.txt
        |-- exonLevelSimilarities-1.0.txt
        |-- org1-ucsc-1.0.bed
        |-- org2-ucsc-1.0.bed
        |-- transcriptLevelMappings-1.0.txt
        `-- transcriptLevelSimilarities-1.0.txt
```

##### (Note: The _exonLevelMappings-1.0.txt_ & _transcriptLevelMappings-1.0.txt_ file contains the mapped exon and transcript pairs from ENSG00000141510-ENSMUSG00000059552 orthologous gene-pair)

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
