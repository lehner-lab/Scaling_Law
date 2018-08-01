# Test all sQTL candidates

In this document, I explain the code found in [006a\_test\_all\_sqtl\_candidates.sh](006a_test_all_sqtl_candidates.sh), where I go through each of the exon skipping events in `Gene_Regions.bed` and perform log-likelihood tests to check whether any snp inside the same gene as the exon of interest may be an sQTL for that splicing event.

## 1. Slow code

Before anything else, we'll store the absolute locations of two of the files produced in earlier steps:

```bash
# store the absolute location of two files we'll use in this script
GENE_REGIONS_BED=absolute_path_to_Gene_Regions.bed_file_generated_earlier
PRUNED_GENOTYPES_FILE=absolute_path_to_0.2.Pruned.genotypes.MAF01.vcf.gz_file_generated_earlier
```
We'll also create a folder called `sQTL_tests` where we'll save all the output from this section:

```bash
# move to directory where I'll keep all the information
mkdir sQTL_tests
cd sQTL_tests
```

There are 35153 exon skipping events in `Gene_Regions.bed` that we can study, so we'll run a for loop where the code will iterate through each of these alternative splicing events, and perform a statistical test for the appropriate sQTL candidates:

```bash
# iterate through the 35153 exon skipping events
for this_number in `seq 1 35153`;
do
    # code here
done
```
The code that goes inside the for loop is as follows. First, we'll create a unique folder for this particular exon skipping event:

```bash
# make a new directory and move to it
mkdir $this_number
cd $this_number
```
Subset `Gene_Regions.bed` to only include the row corresponding to the splicing event of interest:

```bash
# subset Gene_Regions.bed file to only include the gene we're interested in right now
tail -n +$(( this_number + 1 )) $GENE_REGIONS_BED  | head -n 1 > ./Gene_Regions.bed
```
Similarly, filter the VCF file so that it only includes variants inside the gene where the exon skipping event is located:

```bash
# now subset the vcf file and extract it
bcftools view --regions-file ./Gene_Regions.bed $PRUNED_GENOTYPES_FILE -O z -o ./genotypes.MAF01.GeneRegions.vcf.gz
gunzip ./genotypes.MAF01.GeneRegions.vcf.gz
```
The filtered VCF file is then parsed with [006b\_process\_vcf.awk](006b_process_vcf.awk), which returns a table with the following format:

| SNP/Variant | Sample 1 | Sample 2 | Sample 3 | ... |
|-------------|----------|----------|----------|-----|
| SNP ID 1    | 0        | 1        | 1        | ... |
| SNP ID 2    | 1        | 0        | 0        | ... |
| SNP ID 3    | 0        | 1        | 2        | ... |
| ...         | ...      | ...      | ...      | ... |

where 0 means the snp/variant is not found in the sample, 1 means there is one copy of the snp present in the sample DNA (i.e. this sample is heterozygous for this snp), 2 means there are two copies of the snp present (i.e. homozygous for the snp).

```bash
# process the vcf file
awk -f ../../006b_process_vcf.awk ./genotypes.MAF01.GeneRegions.vcf > ./genotypes.MAF01.GeneRegions.Parsed.txt
```
Now, run [006c\_sQTL\_test.R](006c_sQTL_test.R) to perform likelihood-ratio test for all the variants that could be a potential sQTL for this particular alternative splicing event:

```bash
# calculate p values for all possible qtls
Rscript ../../006c_sQTL_test.R 
```
Finally move back up before reating the next folder:

```bash
# move back up
cd ..
```

## 2. Faster code



Alternatively, if you have access to a cluster with an SGE batch system, you can simply submit the following code to the queue. A separate job will be sent for each of the 35153 exon skipping events, which means the whole thing will take far less time:

```bash
#!/bin/bash
#$ -t 1-35153
#$ -o location_of_queue_output_files
#$ -e location_of_queue_error_files

# save SGE_TASK_ID variable with more user-friendly name
this_number=${SGE_TASK_ID}

# store the absolute location of two files we'll use in this script
GENE_REGIONS_BED=absolute_path_to_Gene_Regions.bed_file_generated_earlier
PRUNED_GENOTYPES_FILE=absolute_path_to_0.2.Pruned.genotypes.MAF01.vcf.gz_file_generated_earlier

# move to sQTL_tests  directory where I'll keep all the information
# (must have been created before submitting this job to the queue)
cd absolute_path_to_sQTL_tests

# make a new directory and move to it
mkdir $this_number
cd $this_number

# subset Gene_Regions.bed file to only include the gene we're interested in right now
tail -n +$(( this_number + 1 )) $GENE_REGIONS_BED  | head -n 1 > ./Gene_Regions.bed

# now subset the vcf file and extract it
bcftools view --regions-file ./Gene_Regions.bed $PRUNED_GENOTYPES_FILE -O z -o ./genotypes.MAF01.GeneRegions.vcf.gz
gunzip ./genotypes.MAF01.GeneRegions.vcf.gz

# process the vcf file
awk -f ../../006b_process_vcf.awk ./genotypes.MAF01.GeneRegions.vcf > ./genotypes.MAF01.GeneRegions.Parsed.txt

# calculate p values for all possible qtls
Rscript ../../006c_sQTL_test.R 
```