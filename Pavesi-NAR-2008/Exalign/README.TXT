Exalign 1.0

Exalign is a tool that performs alignments of gene exonic structures. 
In particular, it works by aligning arrays representing the exon lengths of genes. 
While the aligning algorithms used by Exalign are the well known Smith-Waterman and Needleman-Wunsch ones (with a suitable scoring function), the program 
performs a post-processing step that tries to identify intron gain/loss events. For further info on how the alignment algorithm works please refer to the article.

In this file you find help and examples on how to use the standalone commandline version of Exalign.

Please note that  Exalign needs at least one file containing genes structure informations. The file must have the same format of the ones you should have found in this package: Homo_sapiens.rf and Mus_musculus.rf.
The tab separated fields of each line are: transcript id, chromosome, strand, genomic cds start, genomic cds end, number of exons, genomic position of all exon starts (comma separated), genomic position of all exon ends (comma separated), gene name (optional), protein product id (optional).  
You can easily download files formatted in this way by using the UCSC table browser at http://genome.ucsc.edu/cgi-bin/hgTables for every organism available in the database, or by using the Exalign page from which you downloaded this package.

TO COMPILE:

Enter the exalign directory, and type "./compile". Exalign is written in C++, and if an error message appears, it's probably because "gcc/g++" is not installed on your computer. Edit the "compile" file, and replace "g++" with the name of the C++ compiler you have. If no error message appears, then the executable called "exalign" has been produced. You're free to move it wherever you want, but in order to access genome annotations/frequency files, you have to keep all the "rf" and "freq" files in the same directory (folder) of the executable. 

Exalign can be run in different modes:

a given gene (query) is aligned against another gene (target)
a given gene (query) is compared against a collection of genes (target)
all genes of a collection (query) are compared against all genes of another collection (target)

Each comparison (alignment) can be performed in three modes: local, global, or "glocal" (see the article). 

To run the program, type "./exalign" followed by one or more of the following options:

ALIGNMENT MODE: default is "global". Add "-l" for local alignment, or "-gl" for "glocal" alignment

QUERY: 

You must specify a file containing query gene structures:

-q file_name: The "query file". This file must contain one or more gene structure(s). Exalign will align every gene structure contained in the "query file" against every gene structure contained in the "organism file" (unless the -qs and/or -gs options are specified). For every gene in the query file will be reported only the best alignments obtained against the genes contained in the "organism file" (see below).

-qs gene_id: Using this option you choose to align only the selected "gene_id" gene structure from your "query file".

TARGET:

You must specify the file containing the "target" gene structure(s).

-O file_name: The "organism file". This file must contain one or more gene structure(s). You can look at this file as the database you are querying using the gene(s) contained in your query file.  Please note that the "organism file" and the "query file" can be the same file.

-gs gene_id: Using this option you choose to use only the selected "gene_id" (target) gene structure from your "organism file".

OTHER OPTIONS:


-freq This option will calculate the frequencies for the "organism file" in use. You need to use it once only the first time you use a new "organism file". Please note that it's warmly recommended to use, as an "organism file", a file containing the largest possible number of gene structures for the same organism. This allows for accurate exon length frequency estimation and will assure the best results from Exalign algorithm.

-t A file containing a summary table of the results is written. This option is useful when doing large scale analysis. The file will have the same name of the "query file" with the additional extension .tab. 

-og num: Set the OPENING GAP penalty to num. Since Exalign uses a dynamic gap penalty procedure you do not really need to change this value. Changing this value will moreover affect the statistical significance of the E-values computed by Exalign. 

-cg num: Set the CONTINUING GAP penalty to num. Since Exalign uses a dynamic gap penalty procedure you do not really need to change this value. Changing this value will moreover affect the statistical significance of the E-values computed by Exalign.

-ss num: Set the number of start exons that are not aligned. The default value of one is a good choice in almost every case.

-es num: Set the number of end exons that are not aligned. The default value of one is a good choice in almost every case.

-el num: Set how many aligned exons are displayed per line in the output. Change this value only if your output is messed up.

-h displays a quick help.

DATABASE SEARCH SPECIFIC OPTIONS:

-b num: How many results you want to be displayed in your output whenever you align a gene against many genes (database search). Please note that if that value exceeds the one set with -B it will be set to the same value of -B.

-B num: How many of the highest scoring alignments will be subjected to the merging exons step (the step that allows intron gain/loss detection). Increasing this number usually permits to investigate for intron gain/loss events in distant homologues but it may greatly increase the required computational time.

EXAMPLES:

*******  DATABASE SEARCH (1)

Suppose you want to search the Human gene structure database for genes with a structure similar to the one of the TP53 refseq transcript using the default global alignment algorithm.

Then you can type:

./exalign -q Homo_sapiens.rf -qs NM_000546 -O Homo_sapiens.rf 

That is - from the file containing all human RefSeq annotations, find gene NM_000546, and compare it against all genes of the same set.

The first result you obtain will look like that:

Query = NM_000546 Gene = NM_000546


(223)   102     22      279     184     113     110     137     74      107     (1289)  (-)
        |       |       |       |       |       |       |       |       |
(223)   102     22      279     184     113     110     137     74      107     (1289)  (-)

Exact matches:          9 on 9 (100%)
Nearly exact matches:   0 on 9 (0%)
Other matches:          0 on 9 (0%)

Score = 44.39 ( from 44.39 Avg.Score = 0 Dev.Std = 0 )
Evalue = 1.335e-15
**************


The top array represents the exon lengths of the query gene, while the bottom array represents the exon lengths of the highest scoring gene structure. Naturally, since we searched the Human database using a human gene, the highest scoring structure coincides with the query gene. The first and the last exon of both gene structures are enclosed in parentheses to underline the fact that they were not actually aligned and are output only for completeness. The '|' symbols represent exact matches betweeen exon lenghts.

The second result you obtain will look like this:

Query = NM_000546 Gene = NM_005427


(223)   102     22      279     184     113     110     137     74      107     -
        :.      *.      *.      :.      :.      |       :.      *.      *.
(77)    98      121     243     187     116     110     143     89      122     149

        -       -       (1289)  (-)

        139     94      (546)   (+)

Exact matches:          1 on 9 (11.11%)
Nearly exact matches:   4 on 9 (44.44%)
Other matches:          4 on 9 (44.44%)
Score = 17.28 ( from 17.28 Avg.Score = 0 Dev.Std = 0 )
Evalue = 0.0007888
**************


The second most similar gene structure to the one of NM_000546 found in the Human database is NM_0005427, that is TP73. The symbol ":" in the alignment represents a nearly exact match (length difference <= 12), while the "*" represents a mismatch. When the "." appears near those two symbols it means that the two aligned exons share the same frame (see the article for more info about frames).

Query = NM_000546 Gene = NM_003722


(223)   102     22      279     184     113     110     137     74      107     -
        *       *.      *.      :.      :.      |       |       :.      *.
(151)   129     133     255     187     116     110     137     83      137     158

        -       -       (1289)  (-)

        145     94      (3068)  (+)

Exact matches:          2 on 9 (22.22%)
Nearly exact matches:   3 on 9 (33.33%)
Other matches:          4 on 9 (44.44%)

Score = 16.04 ( from 16.04 Avg.Score = 0 Dev.Std = 0 )
Evalue = 0.002749

The third most similar gene structure to the one of NM_000546 found in the Human database is NM_003722, or TP73L. When an exon is aligned with a "-" it means that there is a gap in the alignment.

There will be two more alignments in your output since the default value for the -b option is 5.

******* EXAMPLE OF LOCAL ALIGNMENT 

When you choose the local alignment algorithm, instead of the default global one, there will be a number between brackets near every exon length. This number represents the position of the exon in the gene structure.

Suppose you want to align the gene structure of NM_000546 (Human TP53) with NM_011641 (Mouse Trp63) using local alignment algorithm.
You should type:

./exalign -q Homo_sapiens.rf -qs NM_000546 -O Mus_musculus.rf -gs NM_011641 -l

and the output will be:

Query = NM_000546 Gene = NM_011641


279[4]  184[5]  113[6]  110[7]  137[8]  74[9]   107[10]
*.      :.      :.      |       |       :.      *.
255[2]  187[3]  116[4]  110[5]  137[6]  83[7]   137[8]

Exact matches:          2 on 7 (28.57%)
Nearly exact matches:   3 on 7 (42.86%)
Other matches:          2 on 7 (28.57%)

Score = 21.13 ( from 21.13 )
Evalue = 0.0009923
**************

Note that 279 is the length of the fourth exon of NM_000546, 184 of the fifth and so on.

******* 1 VS 1 ALIGNMENT

When a M appears near an exon it means that this exon is the result of the merging routine thought to detect intron gain/loss events . 

For example, let align the Human TP53 with NM_030989 (Rat Tp53).

./exalign -q Homo_sapiens.rf -qs NM_000546 -O Rattus_norvegicus.rf -gs NM_030989

The output will be something like this:

Query = NM_000546 5(113)+6(110) Gene = NM_030989


(223)   102     22      279     184     223M    137     74      107     (1289)  (-)
        *.      |       :.      |       |       |       |       |
(175)   83      22      273     184     223     137     74      107     (507)   (+)

Exact matches:          6 on 8 (75%)
Nearly exact matches:   1 on 8 (12.5%)
Other matches:          1 on 8 (12.5%)

Score = 34.7 ( from 28.37 )
Evalue = 1.212e-11
**************

The M near the 223 in the fifth aligned exon of the query gene means that this is a result of the merging step. Note that the reported name of the query gene is now changed in NM_000546 5(113)+6(110) to reflect this fact. That name means that the 5th (113bp long) and the 6th (110bp long) aligned exons of the original alignment have been merged to produce a better alignment.  Note also that the score of the alignment is now 34.7, while 28.37 was the score of the alignment before the union of the two exons.  When an M appears in your alignment it could be worth to check with more attention if what you are looking at is an intron gain/loss event.

******* USING NEW STRUCTURE FILES

Suppose you want to work with gene structures of your favourite species and you downloaded them from UCSC genome browser in a file called Bos_taurus.rf. What you have to do before using it is to type:

exalign -O Bos_taurus.rf -freq

The file called Bos_taurus.rf.freq will be generated and from now on you can work with Bos_taurus.rf.

******* LARGE SCALE ANALYSIS

Now you want to align all the gene structures contained in Bos_taurus.rf with all the Human ones to perform a genome wide analysis. You want that only the best result for every alignment is reported too.

Type:

./exalign -q Bos_taurus.rf -O Homo_sapiens.rf -t -b 1

Now exalign will align every gene in Bos_taurus.rf with every gene in Homo_sapiens.rf reporting only the highest scoring alignment for each gene of Bos_taurus.rf (it may take a lot of time!). A file called Bos_taurus.rf.tab will be also generated, containing a summary of your results. 

CONTACTS

For every trouble you may encounter using Exalign please feel free to contact us at:

federico.zambelli@unimi.it or giulio.pavesi@unimi.it


