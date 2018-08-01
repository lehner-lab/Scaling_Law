# Identify common subjects

Not all GTEx subjects were genotyped and not all had their splicing junction reads counted. Therefore, to carry out an analysis involving both junctions and genotypes, I need to find out which subjects are common to both datasets. In this document, I explain the code found in [002\_common\_subjects.sh](002_common_subjects.sh), which is written in bash.


## 1. Download junctions file

I downloaded the junction read counts file (`GTEx_Analysis_2016-01-15_v7_STARv2.4.2a_junctions.gct.gz`) from the GTEx Project [downloads page](https://www.gtexportal.org/gtex_analysis_v7/datasets), which is freely accessible to any user who has a GTEx login account.


## 2. Subjects in both the junctions and genotype files

We need to take all the subjects present in the genotype file, and all the subjects present in the junction read counts file, and find the intersect. We first extract the subjects found in the read counts file:

```bash
# get subjects from junction read counts file
gunzip GTEx_Analysis_2016-01-15_v7_STARv2.4.2a_junctions.gct.gz
tail -n +3 GTEx_Analysis_2016-01-15_v7_STARv2.4.2a_junctions.gct | head -n 1 | awk 'BEGIN{OFS="\n"}{$1=$1; print $0}' | tail -n +3 | cut -d'-' -f 2 | awk '{print "GTEX-"$0}' | sort | uniq > Junctions_Subjects.txt
```
Then, the subjects in the genotype file:

```
# get subjects from genotypes file
gunzip genotypes.MAF01.vcf.gz
grep -v "^##" genotypes.MAF01.vcf | head -n 1 | awk 'BEGIN{OFS="\n"}{$1=$1; print $0}' | tail -n +4 | sort | uniq > Genotypes_Subjects.txt
```
Finally, which subjects are found in both files?

```
# intersect of both files
comm -12 Junctions_Subjects.txt Genotypes_Subjects.txt > Common_Subjects.txt
```