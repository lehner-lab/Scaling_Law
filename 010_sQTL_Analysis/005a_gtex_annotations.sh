#!/bin/bash

# pre-process sample annotations file
cut -f1,6,9 GTEx_v7_Annotations_SampleAttributesDS.txt | awk 'BEGIN {FS="\t";OFS="\t"}; {gsub(/-/, ".", $1)}{print $0}' > IDs_Tissues_IschaemicTime.txt

# pre-process subject annotations file
grep -v -e '^#' -e '^$' phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt | head -n -1 | cut -f 2,4,5 | sed s/-/\./ > Subject_Sex_Age.txt