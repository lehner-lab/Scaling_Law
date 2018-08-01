# Download & filter genotype file

In this document, I explain the code in [001\_filter\_genotypes\_file.sh](001_filter_genotypes_file.sh). To analyse the dependence of sQTL effects on the starting level of exon inclusion, we chose to analyse data from the [GTEx Project](https://www.gtexportal.org/home/). This allowed us to look at the effect of the same sQTL in different tissues (i.e. different 'starting PSIs' if the exon studied is differentially included in different tissues).

Please note that many of the intermediate files produced in this document contain sensitive information extracted from a file whose access is restricted. Furthermore, many of these intermediate files are very big. Therefore, I have not uploaded them to this folder. The instructions here should be enough to reproduce the missing files.

The code in this file is written in bash, and to run it you will need:

* **BCFTools v1.6** installed (downloadable from [here](https://sourceforge.net/projects/samtools/files/samtools/1.6/) although the [latest](https://samtools.github.io/bcftools/) version should also work).

## 1. Download file

The GTEx genotype file is protected and you'll need to apply for access through [dbgap](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v7.p2):

* Genotype matrix in VCF format: `phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1.GRU.tar` (the file we're interested in is `GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz` and is located inside this massive tar file)


## 2. Filter for variants common in the population

The genotype calls file is huge. I subsetted it to include only potential QTLs that are common in the population. The first thing I did was to rename the file and move it to the `Data` folder:

```bash
# rename and move to Data/ folder
mv GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz genotypes.vcf.gz
mv genotypes.vcf.gz Data/
```
Next, I subsetted the VCF file to include only potential QTLs that occur in at least 1% of subjects in our cohort. To do this, I used [BCFtools](https://samtools.github.io/bcftools/):

```bash
# subset the file to only include common snps, indels etc
bcftools view --min-af 0.01:minor Data/genotypes.vcf.gz -O z -o Data/genotypes.MAF01.vcf.gz  # 5-6 hours! from 120gb to 37gb
```
**Warning**: depending on how the VCF file was handled, you might get an error further downstream, saying something like *Could not load the index*. This error has been [reported before](https://github.com/samtools/bcftools/issues/129) and it is caused by your file having been compressed with gzip instead of [bgzip](http://www.htslib.org/download/). The fix is as follows:

```bash
# error reported in https://github.com/samtools/bcftools/issues/129
# the fix is this following line
zcat genotypes.MAF01.vcf.gz | bgzip -c > NEW.genotypes.MAF01.vcf.gz && tabix NEW.genotypes.MAF01.vcf.gz # this step can take 3~4 hours

# and rename
rm genotypes.MAF01.vcf.gz
mv NEW.genotypes.MAF01.vcf.gz genotypes.MAF01.vcf.gz
```


## 3. Remove variants in linkage disequilibrium

Many variants may co-occur and be associated with the exact same phenotypic changes. In this case, they are all effectively the same QTL. To remove these redundant variants, I used `bcftools +prune` with an R2 limit of 0.2 within a window of 100kb (default window size).

```bash
# prune VCF file
bcftools +prune --max-LD 0.2 Data/genotypes.MAF01.vcf.gz --output-type z --output Data/Pruned.genotypes.MAF01.vcf.gz 

# rename
rm Data/genotypes.MAF01.vcf.gz
mv Data/Pruned.genotypes.MAF01.vcf.gz Data/genotypes.MAF01.vcf.gz
```