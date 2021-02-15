BALSA: Bayesian algorithm for local sequence alignment

The SmithWaterman algorithm yields a single alignment,
which, albeit optimal, can be strongly affected
by the choice of the scoring matrix and the gap
penalties. Additionally, the scores obtained are
dependent upon the lengths of the aligned
sequences, requiring a post-analysis conversion. To
overcome some of these shortcomings, we developed
a Bayesian algorithm for local sequence alignment
(BALSA), that takes into account the uncertainty
associated with all unknown variables by incorporating
in its forward sums a series of scoring matrices, gap
parameters and all possible alignments. The algorithm
can return both the joint and the marginal optimal
alignments, samples of alignments drawn from the
posterior distribution and the posterior probabilities
of gap penalties and scoring matrices. Furthermore,
it automatically adjusts for variations in sequence
lengths. BALSA was compared with SSEARCH, to
date the best performing dynamic programming
algorithm in the detection of structural neighbors.
Using the SCOP databases PDB40D-B and PDB90D-B,
BALSA detected 19.8 and 41.3% of remote homologs
whereas SSEARCH detected 18.4 and 38% at an error
rate of 1% errors per query over the databases,
respectively

5 directories, 36 files

The bin directory contains BALSAfile and BALSAdatabase compiled with
C++ gcc version 7.3.0 for Ubuntu 20.04.2 LTS.

To test:
bin/BALSAfile -SAMPLES -INFILE1=data/1.fa  -INFILE2=data/2.fa -OUTFILE=out.txt -MATRIX=BLOSUM_30,-10,-5


If used in a publication, please reference
Webb B-JM, Liu JS, Lawrence CE (2002) BALSA: Bayesian algorithm for local
sequence alignment. Nucleic Acids Research 30: 12681277
