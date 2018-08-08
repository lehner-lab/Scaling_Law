#!/bin/bash

# Make SRR-to-GSM table
wget ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab
grep ^SRR SRA_Accessions.tab | grep GSM | awk 'BEGIN {FS="\t"; print "GSM" FS "SRR"}; {print $10 FS $1}' > GSM_SRR.txt
rm SRA_Accessions.tab

# prefetch SRR files from GEO
while IFS=$'\t' read -r -a myArray
do
	SRR_ID=$(grep ${myArray[0]} GSM_SRR.txt | cut -f 2)
	prefetch -v $SRR_ID
done < Data/Sample_IDs.txt

# make a directory to store the fastq files
mkdir fastq_files

# extract fastq files
while IFS=$'\t' read -r -a myArray
do
	SRR_ID=$(grep ${myArray[0]} GSM_SRR.txt | cut -f 2)
	fastq-dump --outdir fastq_files/ --split-files ~/ncbi/public/sra/${SRR_ID}.sra
done < Data/Sample_IDs.txt

# rename files
while IFS=$'\t' read -r -a myArray
do
	SRR_ID=$(grep ${myArray[0]} GSM_SRR.txt | cut -f 2)
	mv fastq_files/${SRR_ID}_1.fastq fastq_files/${myArray[1]}_1.fastq
	mv fastq_files/${SRR_ID}_2.fastq fastq_files/${myArray[1]}_2.fastq
done < Data/Sample_IDs.txt