#!/bin/bash

# Before running this code, you should have downloaded the file named
# GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz
# from dbGap

# rename and move to Data/ folder
mv GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz genotypes.vcf.gz
mv genotypes.vcf.gz Data/

# subset the file to only include common snps, indels etc
bcftools view --min-af 0.01:minor Data/genotypes.vcf.gz -O z -o Data/genotypes.MAF01.vcf.gz  # 5-6 hours! from 120gb to 37gb


# the following code is commented out;
# only necessary if you run into the "Could not load the index" error

# # error reported in https://github.com/samtools/bcftools/issues/129
# # the fix is this following line
# zcat genotypes.MAF01.vcf.gz | bgzip -c > NEW.genotypes.MAF01.vcf.gz && tabix NEW.genotypes.MAF01.vcf.gz # this step can take 3~4 hours
# 
# # and rename
# rm genotypes.MAF01.vcf.gz
# mv NEW.genotypes.MAF01.vcf.gz genotypes.MAF01.vcf.gz


# prune VCF file
bcftools +prune --max-LD 0.2 Data/genotypes.MAF01.vcf.gz --output-type z --output Data/Pruned.genotypes.MAF01.vcf.gz 

# rename
rm Data/genotypes.MAF01.vcf.gz
mv Data/Pruned.genotypes.MAF01.vcf.gz Data/genotypes.MAF01.vcf.gz