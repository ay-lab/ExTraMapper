## Steps to generate the input files
The users should run the _run_human_monkey_preprocess_data_steps.sh_ to generate the inputfiles. All the input files will be generated under _preprocess/data_ folder. The _run_human_monkey_preprocess_data_steps.sh_ has six individual steps and should be run one after the another like the following 
```bash
Set the following path
export EXTRAMAPPER_DIR=/path/to/this/folder

This will fetch organism specific chromosomal fasta, gtf and liftOver files. 
$ ./run_human_mouse_preprocess_data_steps.sh 0


The above step may produce an error while downloading the monkey genome from UCSC. 
In that case, please do the following and the script will produce the required fasta files.
$ cd ./preprocess/data/reference_genomes/rheMac10/
$ perl getFasta.pl
$ cd -

The next step will create the genomedata object files. This step requires genomedata package
which can be installed by running the following commnand and followed by the step 1 commnad.
$ pip install genomedata --user
$ ./run_human_mouse_preprocess_data_steps.sh 1


Step2 will create pickle files and gene summaries
$ ./run_human_mouse_preprocess_data_steps.sh 2


Step3 will run the liftOver with multiple mappings and also compute intersections with the other set of exons
$ ./run_human_mouse_preprocess_data_steps.sh 3


Step 4 to 6 will generate the input files
$ ./run_human_mouse_preprocess_data_steps.sh 4
$ ./run_human_mouse_preprocess_data_steps.sh 5
$ ./run_human_mouse_preprocess_data_steps.sh 6
```
##### The whole process should take some time to finish!
##### (Check also the Human-Mouse data processing steps)
