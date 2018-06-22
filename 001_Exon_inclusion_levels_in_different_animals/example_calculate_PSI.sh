#!/bin/bash 

# download chim heart RNA-seq files
wget http://www.nhprtr.org/data/2014_NHP_tissuespecific/10NHP_Illumina_Baylor/Chimpanzee/Chimpanzee_Heart_359_D2EKWACXX-8-ID22_1_sequence.txt.bz2
wget http://www.nhprtr.org/data/2014_NHP_tissuespecific/10NHP_Illumina_Baylor/Chimpanzee/Chimpanzee_Heart_359_D2EKWACXX-8-ID22_2_sequence.txt.bz2

# extract the files
bunzip2 Chimpanzee_Heart_359_D2EKWACXX-8-ID22_1_sequence.txt.bz2 Chimpanzee_Heart_359_D2EKWACXX-8-ID22_2_sequence.txt.bz2

# build a directory to store STAR output
mkdir STAR_Output

# run STAR
STAR --runThreadN 10  --genomeDir ./genomeDir --readFilesIn Chimpanzee_Heart_359_D2EKWACXX-8-ID22_1_sequence.txt Chimpanzee_Heart_359_D2EKWACXX-8-ID22_2_sequence.txt --outFileNamePrefix STAR_Output/Pan_troglodytes.CHIMP2.1.4.84_Heart_

# convert the output from STAR into a format
# that is similar to the output from TopHat
# oneliner code taken from this protocol:
# http://onlinelibrary.wiley.com/doi/10.1002/0471142905.hg1116s87/abstract
cd STAR_Output
awk 'BEGIN{OFS="\t"}{print $1, $2-20-1, $3+20, "JUNCBJ"NR, $7, ($4 == 1)? "+":"-",$2-20-1, $3+20, "255,0,0", 2, "20,20", "0,300" }' *SJ.out.tab > junctions.bed

# convert STAR .sam output into .bam
samtools view -hb *Aligned.out.sam > accepted_hits.bam

# go back up
mv ..

# make a directory
mkdir PSI_Calculation

# move some files...
mv STAR_Output/accepted_hits.bam PSI_Calculation/
mv STAR_Output/junctions.bed PSI_Calculation/

# move to relevant folder
cd PSI_Calculation

# start counting
bash ../PSI.sh StartPSI ../Exonic_Parts.gtf 101 accepted_hits.bam junctions.bed Chimp_Heart
