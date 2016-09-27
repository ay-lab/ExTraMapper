This entry describes how we downloaded and processed data from (Exalign):

Pavesi, G., Zambelli, F., Caggese, C. and Pesole, G. (2008) Exalign: a new
method for comparative analysis of exon-intron gene structures. Nucleic acids
research, 36, e47.

Through the five steps (0-4) described in runall script we were able to:

- step0: Download source code and compile it
- step1: Convert our input files for ExTraMapper to appropriate format for
  Exalign
- step2: Split the big input file into a file pair per gene pair
- step3: Run exalign on each gene pair (similar to how ExTraMapper runs)
- step4: Combine results of transcript mappings from per gene pair runs

You will need to run all five steps sequentially and some might be 
time/resource consuming (especially step 3):
./runall 0
./runall 1
./runall 2
./runall 3
./runall 4

