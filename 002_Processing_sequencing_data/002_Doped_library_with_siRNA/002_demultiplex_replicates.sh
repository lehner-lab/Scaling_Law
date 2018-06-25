#!/bin/bash

# NOTE: you don't need to run this step if you already downloaded the
# demultiplexed files.

# make and move to folder where output will be deposited
mkdir 002_demultiplexed
cd 002_demultiplexed

# run sabre
# I'm making the following assumptions:
# 1) the raw files are located in the folder just above
#    'demultiplexed' and that they are called Doped_siRNA_1.fastq
#    and Doped_siRNA_2.fastq
# 2) '002_barcodes.txt' is located in the folder above 'demultiplexed'
sabre pe -m 2 -c -f ../Doped_siRNA_1.fastq -r ../Doped_siRNA_2.fastq -b ../002_barcodes.txt -u forward_reads_no_barcode.txt -w reverse_reads_no_barcode.txt

# go back up
cd ..