## Steps to generate the input files
The users should run the _run_human_mouse_preprocess_data_steps.sh_ to generate the inputfiles. All the input files will be generated under _preprocess/data_ folder. The _run_human_mouse_preprocess_data_steps.sh_ has six individual steps and should be run one after the another like the following 
```bash
This will fetch organism specific chromosomal fasta, gtf and liftOver files. 
$ ./run_human_mouse_preprocess_data_steps.sh 0

The next step will create the genomedata object files. This step requires genomedata package
which can be installed by running the following commnand and followed by the step 1 commnad.
$ pip install genomedata --user
$ ./run_human_mouse_preprocess_data_steps.sh 1

