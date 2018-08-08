#!/bin/bash

# Before running this code, you should have downloaded the file named
# GTEx_Analysis_2016-01-15_v7_STARv2.4.2a_junctions.gct.gz from GTEx

# get subjects from junction read counts file
gunzip GTEx_Analysis_2016-01-15_v7_STARv2.4.2a_junctions.gct.gz
tail -n +3 GTEx_Analysis_2016-01-15_v7_STARv2.4.2a_junctions.gct | head -n 1 | awk 'BEGIN{OFS="\n"}{$1=$1; print $0}' | tail -n +3 | cut -d'-' -f 2 | awk '{print "GTEX-"$0}' | sort | uniq > Junctions_Subjects.txt

# get subjects from genotypes file
gunzip genotypes.MAF01.vcf.gz
grep -v "^##" genotypes.MAF01.vcf | head -n 1 | awk 'BEGIN{OFS="\n"}{$1=$1; print $0}' | tail -n +4 | sort | uniq > Genotypes_Subjects.txt

# intersect of both files
comm -12 Junctions_Subjects.txt Genotypes_Subjects.txt > Common_Subjects.txt