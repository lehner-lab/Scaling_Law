#!/bin/bash

# prefetch SRR files from GEO
while IFS=$'\t' read -r -a myArray
do
	SRR_ID=${myArray[0]}
	prefetch -v $SRR_ID
done < Data/Sample_IDs.txt

# make a directory to store the fastq files
mkdir fastq_files

# extract fastq files
while IFS=$'\t' read -r -a myArray
do
	SRR_ID=${myArray[0]}
	fastq-dump --outdir fastq_files/ --split-files ~/ncbi/public/sra/${SRR_ID}.sra
done < Data/Sample_IDs.txt

# rename files
while IFS=$'\t' read -r -a myArray
do
	SRR_ID=${myArray[0]}
	mv fastq_files/${SRR_ID}_1.fastq fastq_files/${myArray[1]}.fastq
done < Data/Sample_IDs.txt