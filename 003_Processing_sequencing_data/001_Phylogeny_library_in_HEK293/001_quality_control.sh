#!/bin/bash

# make a directory for fastqc output
mkdir fastq_reports

# command to run FastQC
# - o argument is output folder which must previously exist
# last argument is the file to be processed by fastqc
# (I'm assuming the fastq files have all been renamed to something
# that starts with "Phylogeny_HEK293_")
fastqc -o ./fastqc_reports Phylogeny_HEK293_*
