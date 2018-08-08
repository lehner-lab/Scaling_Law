#!/bin/bash

# location of PEAR software as a variable
PEAR=path_to_pear_executable


# input files with forward paired-end reads
FORWARD_FILE_1=002_demultiplexed/Control_RNAi_Pair_1_Rep_1.txt
FORWARD_FILE_2=002_demultiplexed/Control_RNAi_Pair_1_Rep_2.txt
FORWARD_FILE_3=002_demultiplexed/Control_RNAi_Pair_1_Rep_3.txt
FORWARD_FILE_4=002_demultiplexed/SF3B_RNAi_Pair_1_Rep_1.txt
FORWARD_FILE_5=002_demultiplexed/SF3B_RNAi_Pair_1_Rep_2.txt
FORWARD_FILE_6=002_demultiplexed/SF3B_RNAi_Pair_1_Rep_3.txt


# input files with reverse paired-end reads
REVERSE_FILE_1=002_demultiplexed/Control_RNAi_Pair_2_Rep_1.txt
REVERSE_FILE_2=002_demultiplexed/Control_RNAi_Pair_2_Rep_2.txt
REVERSE_FILE_3=002_demultiplexed/Control_RNAi_Pair_2_Rep_3.txt
REVERSE_FILE_4=002_demultiplexed/SF3B_RNAi_Pair_2_Rep_1.txt
REVERSE_FILE_5=002_demultiplexed/SF3B_RNAi_Pair_2_Rep_2.txt
REVERSE_FILE_6=002_demultiplexed/SF3B_RNAi_Pair_2_Rep_3.txt


# output files
mkdir 003_merged
MERGED_FILE_1=003_merged/Control_RNAi_Rep_1_Merged
MERGED_FILE_2=003_merged/Control_RNAi_Rep_2_Merged
MERGED_FILE_3=003_merged/Control_RNAi_Rep_3_Merged
MERGED_FILE_4=003_merged/SF3B_RNAi_Rep_1_Merged
MERGED_FILE_5=003_merged/SF3B_RNAi_Rep_2_Merged
MERGED_FILE_6=003_merged/SF3B_RNAi_Rep_3_Merged


# run PEAR software to merge forward and reverse reads into one file
time $PEAR -f $FORWARD_FILE_1 -r $REVERSE_FILE_1 -m 119 -n 119 -v 115 -j 4 -o $MERGED_FILE_1
time $PEAR -f $FORWARD_FILE_2 -r $REVERSE_FILE_2 -m 119 -n 119 -v 115 -j 4 -o $MERGED_FILE_2
time $PEAR -f $FORWARD_FILE_3 -r $REVERSE_FILE_3 -m 119 -n 119 -v 115 -j 4 -o $MERGED_FILE_3
time $PEAR -f $FORWARD_FILE_4 -r $REVERSE_FILE_4 -m 119 -n 119 -v 115 -j 4 -o $MERGED_FILE_4
time $PEAR -f $FORWARD_FILE_5 -r $REVERSE_FILE_5 -m 119 -n 119 -v 115 -j 4 -o $MERGED_FILE_5
time $PEAR -f $FORWARD_FILE_6 -r $REVERSE_FILE_6 -m 119 -n 119 -v 115 -j 4 -o $MERGED_FILE_6


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
