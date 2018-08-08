#!/bin/bash

# location of PEAR software as a variable
PEAR=path_to_pear_executable


# input files with forward paired-end reads
FORWARD_FILE_1=002_demultiplexed/Phylogeny_HEK293_TRBR_Pair_1_Rep_1.txt
FORWARD_FILE_2=002_demultiplexed/Phylogeny_HEK293_TRBR_Pair_1_Rep_2.txt
FORWARD_FILE_3=002_demultiplexed/Phylogeny_HEK293_TRBR_Pair_1_Rep_3.txt
FORWARD_FILE_7=002_demultiplexed/Phylogeny_HEK293_BR_Pair_1_Rep_4.txt
FORWARD_FILE_8=002_demultiplexed/Phylogeny_HEK293_BR_Pair_1_Rep_5.txt
FORWARD_FILE_9=002_demultiplexed/Phylogeny_HEK293_BR_Pair_1_Rep_6.txt
FORWARD_FILE_10=002_demultiplexed/Phylogeny_HEK293_BR_Pair_1_Rep_7.txt
FORWARD_FILE_11=002_demultiplexed/Phylogeny_HEK293_BR_Pair_1_Rep_8.txt
FORWARD_FILE_12=002_demultiplexed/Phylogeny_HEK293_BR_Pair_1_Rep_9.txt


# input files with reverse paired-end reads
REVERSE_FILE_1=002_demultiplexed/Phylogeny_HEK293_TRBR_Pair_2_Rep_1.txt
REVERSE_FILE_2=002_demultiplexed/Phylogeny_HEK293_TRBR_Pair_2_Rep_2.txt
REVERSE_FILE_3=002_demultiplexed/Phylogeny_HEK293_TRBR_Pair_2_Rep_3.txt
REVERSE_FILE_7=002_demultiplexed/Phylogeny_HEK293_BR_Pair_2_Rep_4.txt
REVERSE_FILE_8=002_demultiplexed/Phylogeny_HEK293_BR_Pair_2_Rep_5.txt
REVERSE_FILE_9=002_demultiplexed/Phylogeny_HEK293_BR_Pair_2_Rep_6.txt
REVERSE_FILE_10=002_demultiplexed/Phylogeny_HEK293_BR_Pair_2_Rep_7.txt
REVERSE_FILE_11=002_demultiplexed/Phylogeny_HEK293_BR_Pair_2_Rep_8.txt
REVERSE_FILE_12=002_demultiplexed/Phylogeny_HEK293_BR_Pair_2_Rep_9.txt


# output files
mkdir 003_merged
MERGED_FILE_1=003_merged/Phylogeny_HEK293_TR_Rep_1_Merged
MERGED_FILE_2=003_merged/Phylogeny_HEK293_TR_Rep_2_Merged
MERGED_FILE_3=003_merged/Phylogeny_HEK293_TR_Rep_3_Merged
MERGED_FILE_4=003_merged/Phylogeny_HEK293_BR_Rep_1_Merged
MERGED_FILE_5=003_merged/Phylogeny_HEK293_BR_Rep_2_Merged
MERGED_FILE_6=003_merged/Phylogeny_HEK293_BR_Rep_3_Merged
MERGED_FILE_7=003_merged/Phylogeny_HEK293_BR_Rep_4_Merged
MERGED_FILE_8=003_merged/Phylogeny_HEK293_BR_Rep_5_Merged
MERGED_FILE_9=003_merged/Phylogeny_HEK293_BR_Rep_6_Merged
MERGED_FILE_10=003_merged/Phylogeny_HEK293_BR_Rep_7_Merged
MERGED_FILE_11=003_merged/Phylogeny_HEK293_BR_Rep_8_Merged
MERGED_FILE_12=003_merged/Phylogeny_HEK293_BR_Rep_9_Merged


# run PEAR software to merge forward and reverse reads into one file
$PEAR -f $FORWARD_FILE_1 -r $REVERSE_FILE_1 -m 116 -n 116 -v 116 -j 4 -o $MERGED_FILE_1
$PEAR -f $FORWARD_FILE_2 -r $REVERSE_FILE_2 -m 116 -n 116 -v 116 -j 4 -o $MERGED_FILE_2
$PEAR -f $FORWARD_FILE_3 -r $REVERSE_FILE_3 -m 116 -n 116 -v 116 -j 4 -o $MERGED_FILE_3
$PEAR -f $FORWARD_FILE_1 -r $REVERSE_FILE_1 -m 119 -n 119 -v 115 -j 4 -o $MERGED_FILE_4
$PEAR -f $FORWARD_FILE_2 -r $REVERSE_FILE_2 -m 119 -n 119 -v 115 -j 4 -o $MERGED_FILE_5
$PEAR -f $FORWARD_FILE_3 -r $REVERSE_FILE_3 -m 119 -n 119 -v 115 -j 4 -o $MERGED_FILE_6
$PEAR -f $FORWARD_FILE_7 -r $REVERSE_FILE_7 -m 119 -n 119 -v 115 -j 4 -o $MERGED_FILE_7
$PEAR -f $FORWARD_FILE_8 -r $REVERSE_FILE_8 -m 119 -n 119 -v 115 -j 4 -o $MERGED_FILE_8
$PEAR -f $FORWARD_FILE_9 -r $REVERSE_FILE_9 -m 119 -n 119 -v 115 -j 4 -o $MERGED_FILE_9
$PEAR -f $FORWARD_FILE_10 -r $REVERSE_FILE_10 -m 119 -n 119 -v 115 -j 4 -o $MERGED_FILE_10
$PEAR -f $FORWARD_FILE_11 -r $REVERSE_FILE_11 -m 119 -n 119 -v 115 -j 4 -o $MERGED_FILE_11
$PEAR -f $FORWARD_FILE_12 -r $REVERSE_FILE_12 -m 119 -n 119 -v 115 -j 4 -o $MERGED_FILE_12


####################################
# ABOUT THE ARGUMENTS USED BY PEAR #
####################################

# -f Name of the file containing the forward paired-end reads
# -r File containing reverse-end reads
# -o Output file
# -v Minimum overlap size
# -m Max length of the assembled sequences
# -n Min length of the assembled sequences
# -j Number of threads to use (run in parallel)
# -p P-value for the statistical test (-p 1 disables the test, and this is useful in case of a small overlap; in the case of a large overlap with good sequence quality, it can be set to 0.05)
