#!/bin/bash


# directory where I'll save output from this step
mkdir 004_counts



####################################
#### Software Used in this Step ####
####################################


# fastq-grep (part of fastq-tools): given a regex, print all FASTQ entries
# with a matching nucleotide sequence.
FASTQ_GREP=path_to_fastq-grep

# fastx_reverse_complement (part of the FASTX-Toolkit): produce the RC of
# each entry in a FASTA/FASTQ file.
FASTX_RC=path_to_fastx_reverse_complement

# fastx_collapser (part of the FASTX-Toolkit): collapses identical sequences
# in a FASTA/Q file into a single sequence, maintaining read counts
FASTX_COLLAPSE=path_to_fastx_collapser

# seqtk: a tool for processing sequences in FASTA/Q format. Can trim nucleotides
# from both ends of the sequence.
SEQTK_TRIM=path_to_seqtk

# custom python script to count mutations
COUNT_VARIANTS=004_count_mutations.py



#####################################
#### Sequences used in this expt ####
#####################################

# reference sequence (human sequence)
WTSEQ=GATCCAGATCTAACTTGGGGTGGCTTTGTCTTCTTCTTTTGCCAATTCCACTAATTGTTTGGG

# primers from the forward sequence
PRIMER=CAGCAACACCAAGTGCAAAGAGGAAG
PRIMER_END=TGAAGAGAAAGGAAGTACAGAAAACATGCA

# primers from the reverse sequence
RC_PRIMER=TGCATGTTTTCTGTACTTCCTTTCTCTTCA
RC_PRIMER_END=CTTCCTCTTTGCACTTGGTGTTGCTG



####################
#### Count!!!!! ####
####################

for Sample in Control SF3B;
do
	for i in `seq 1 3`;
	do
		# file used as input for this step
		INPUT_FILE=003_merged/${Sample}_RNAi_Rep_${i}_Merged.assembled.fastq
		
		# files created in this step 
		FILE_WITH_FORWARD_READS=004_counts/${Sample}_RNAi_Rep_${i}.forward
		FILE_WITH_REVERSE_READS=004_counts/${Sample}_RNAi_Rep_${i}.reverse
		FILE_WITH_REVERSE_READS_RC=$FILE_WITH_REVERSE_READS.rc
		CONCATENATED_FILE=004_counts/${Sample}_RNAi_Rep_${i}.concatenated
		COLLAPSED_FILE=004_counts/${Sample}_RNAi_Rep_${i}.collapsed
		TRIMMED_FILE=004_counts/${Sample}_RNAi_Rep_${i}.trimmed
		FINAL_FILE_WITH_COUNTS=004_counts/${Sample}_RNAi_Rep_${i}.counts
		
		# separate forward and reverse reads
		# use fast-grep to use primers to determine which a forward reads
		# and which are reverse reads and drop them in different files
		# usage: fastq-grep regex INPUT_FILE [STDOUT]
		echo $(date) "Separating forward and reverse reads"
		$FASTQ_GREP ^$PRIMER.+$PRIMER_END$ $INPUT_FILE > $FILE_WITH_FORWARD_READS
		$FASTQ_GREP ^$RC_PRIMER.+$RC_PRIMER_END$ $INPUT_FILE > $FILE_WITH_REVERSE_READS
		
		# use fastx_reverse_complement to reverse complement the file with reverse reads
		# usage: fastx_reverse_complement -i INPUT_FILE
		# [-Q33 (because Sanger / Illumina 1.9 encoding)] [STDOUT]
		echo $(date) "Reverse complementing the file with reverse reads"
		$FASTX_RC -i $FILE_WITH_REVERSE_READS -Q33 > $FILE_WITH_REVERSE_READS_RC
		
		# concatenate file with forward reads with reverse-complemented reverse reads
		echo $(date) "Concatenating file with forward reads and file with reverse-complemented reverse reads"
		cat $FILE_WITH_FORWARD_READS $FILE_WITH_REVERSE_READS_RC > $CONCATENATED_FILE
		
		# use SeqTK to remove primer sequences from the start and end of the sequences
		# usage: seqtk trimfq -b number_of_nucleotides_to_be_removed_BEFORE_sequence -e number_of_nucleotides_to_be_removed_at_the_END INPUT_FILE [STDOUT]
		echo $(date) "Trimming primers"
		$SEQTK_TRIM trimfq -b 26 -e 30 $CONCATENATED_FILE > $TRIMMED_FILE
		
		# use fastx_collapse to collapse concatenated file and remove repeated sequences
		# usage the same as fastx_reverse_complement
		echo $(date) "Collapse concatenated file and remove repeated sequences. Keeping number of occurrences."
		$FASTX_COLLAPSE -i $TRIMMED_FILE -Q33 > $COLLAPSED_FILE
		
		# use python script to count mutations and return a table
		# usage: python count_mutations.py INPUT_FASTA_FILE WT_SEQUENCE [STDOUT]
		echo $(date) "Counting variants... almost there..."
		$COUNT_VARIANTS $COLLAPSED_FILE $WTSEQ > $FINAL_FILE_WITH_COUNTS
		
		echo $(date) "Done."
		
	done
done










