# ExTraMapper
ExTraMapper is a tool to find Exon and Transcript-level Mappings of a given pair of orthologous genes between two organisms using sequence conservation. The figure below shows the overall schematic description of ExTraMapper mapping the homologous transcript and exon-pairs between human and mouse genome. 


![ExTraMapper_Figure](https://user-images.githubusercontent.com/18036388/90572310-8b693e00-e168-11ea-9fbc-8188c2834de9.jpg)

# Steps to run ExtraMapper (For python version 3 or later usage)
### Step 1: Prepare the input files
ExTraMapper requires a set of preprocessed files to find the conservation scores. Examples to create these files are provided within [__Human-Mouse-Preprocess-Data__](https://github.com/ay-lab/ExTraMapper/tree/master/Human-Mouse-Preprocess-Data) and [__Human-Mokey-Preprocessed-Data__](https://github.com/ay-lab/ExTraMapper/tree/master/Human-Monkey-Processed-Data) folders. Users should look into these folders and follow the instructions to create the required input files before going to the next step.   


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
