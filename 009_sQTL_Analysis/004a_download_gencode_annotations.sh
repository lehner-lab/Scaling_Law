#!/bin/bash

# download file
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz

# unpack the file
gunzip gencode.v19.annotation.gff3.gz