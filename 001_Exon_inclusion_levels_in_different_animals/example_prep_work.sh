#!/bin/bash 

# downloade chimp genome sequence
wget ftp://ftp.ensembl.org/pub/release-84/fasta/pan_troglodytes/dna/Pan_troglodytes.CHIMP2.1.4.dna.toplevel.fa.gz

# extract file
gunzip Pan_troglodytes.CHIMP2.1.4.dna.toplevel.fa.gz

# download chimp genome annotations
wget ftp://ftp.ensembl.org/pub/release-84/gtf/pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.84.gtf.gz

# extract the file
gunzip Pan_troglodytes.CHIMP2.1.4.84.gtf.gz

# download chimp genome annotations
wget ftp://ftp.ensembl.org/pub/release-84/gtf/pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.84.gtf.gz

# extract the file
gunzip Pan_troglodytes.CHIMP2.1.4.84.gtf.gz

# build folder to dump STAR indices
mkdir genomeDir

# build indices
STAR --runThreadN 5 --limitGenomeGenerateRAM 80000000000 --runMode genomeGenerate --genomeDir ./genomeDir --genomeFastaFiles Pan_troglodytes.CHIMP2.1.4.dna.toplevel.fa --sjdbGTFfile Pan_troglodytes.CHIMP2.1.4.84.gtf
