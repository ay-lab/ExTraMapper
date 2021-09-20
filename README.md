# ExTraMapper
ExTraMapper is a tool to find Exon and Transcript-level Mappings of a given pair of orthologous genes between two organisms using sequence conservation. The figure below shows the overall schematic description of ExTraMapper mapping the homologous transcript and exon-pairs between human and mouse genome. 


![ExTraMapper_Figure](https://user-images.githubusercontent.com/18036388/90572310-8b693e00-e168-11ea-9fbc-8188c2834de9.jpg)

# Steps to run ExtraMapper (For python version 3 or later usage)

### Step 1: Prepare the input files
ExTraMapper requires a set of preprocessed files to find the conservation scores. Examples to create these files are provided within the following folders
 
1. [__Human-Mouse-Preprocessed-Data__](https://github.com/ay-lab/ExTraMapper/tree/master/Human-Mouse-Preprocess-Data) 

    and 
    
2. [__Human-Rhesus_macaque-Preprocessed-Data__](https://github.com/ay-lab/ExTraMapper/tree/master/Human-Monkey-Processed-Data) 

### Steps to generate the input files
The users should run the _extMpreprocess_ to generate the inputfiles within the above Preprocessed-Data folders. All the input files will be generated under _preprocess/data_ folder. All the required executables and scripts are provided here. The _extMpreprocess_ has 7 individual steps and should be run in the following manner 
 
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
    
    example: 
    
    ./extMpreprocess config.human-mouse.conf all
    ```
   <br>
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

Note: The __exonLevelMappings-1.0.txt__ & __transcriptLevelMappings-1.0.txt__ file contains the mapped exon and transcript pairs from __ENSG00000141510-ENSMUSG00000059552__ orthologous gene-pair. 

<br>

# OR

### Step 3: Run ExTraMapper for all the gene pairs
```bash
$ python ExTraMapper.py -h
usage: ExTraMapper.py [-h] -m MAPPING -o1 ORG1 -o2 ORG2 -p all
```

<br>

### Summarise the ExTraMapper results ###
Run _extMsummarise_ script to generate a concatenated file will all the results. Run the script in the follwoing manner 
```bash
$ ./extMsummarise help
Type ./extMsummarise <preprocess_folder> <extramapper_folder> <orthologous_genepair_list> <org1name> <org2name> <outputprefix>
preprocess_folder  : Path to the preprocess folder generated by the extMpreproces script
extramapper_folder : Path to the output folder generated by ExTraMapper program
orthologous_genepair_list : A list of orthologous gene-pairs
org1name : org1 name e.g. human
org2name : org2 name e.g. mouse
outputprefix : output file prefix

example : 
./extMsummarise ./preprocess ./output gene-pair.list human mouse extramapper-result
```
<br>


# Prepocessed Results

Check the [Result/Exon-Pairs](https://github.com/ay-lab/ExTraMapper/tree/master/Result/Exon-Pairs) and [Result/Transcript-Pairs](https://github.com/ay-lab/ExTraMapper/tree/master/Result/Transcript-Pairs) to download the precomputed ExTraMapper result for human-mouse and human-rhesus orthologous exon and transcript pairs.

<br>

### Refer the work
[_ExTraMapper: Exon- and Transcript-level mappings for orthologous gene pairs._](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btab393/6278896?redirectedFrom=fulltext)

__Chakraborty A, Ay F, Davuluri RV. ExTraMapper: Exon- and Transcript-level mappings for orthologous gene pairs. Bioinformatics. 2021 May 20:btab393. doi: 10.1093/bioinformatics/btab393. Epub ahead of print. PMID: 34014317.__

The data shown in the above paper was performed using Human & Mouse ENSMBL release 81 with python 2.7 code. 
The current update is with ENSMBL release 102 and python 3 or later version. To see the older code and data please
change the __Branch__ to [__ExTraMapper-python2v__](https://github.com/ay-lab/ExTraMapper/tree/ExTraMapper-python2v) from __master__

### Check the webserver for a nice vizualization 
https://ay-lab-tools.lji.org/extramapper/index.html
