# Test all potential sQTL's

In this document I describe how I processed `genotypes.MAF01.vcf.gz` to test whether each variant is significantly associated with changes in the PSI of a particular exon.

This document contains three sections, corresponding to three scripts that were used to call sQTL's. The first one is a **bash script**, which loops over all 4885 alternative exon events to test which sQTL's might be significant. The bash script calls two other scripts: an **awk script** and an **R script**. In each section I'll describe what each of these three scripts is doing.

## Bash script

This script does the following:

* create a new folder for each splicing event studied
* creates a small VCF file containing only variants inside the same gene as the exon studied
* processes the vcf file with the awk script (see below)
* finds the significance of each potential sQTL using the R script (see below)


```bash
for this_number in `seq 1 4885`; do
	# make a new directory and move to it
	mkdir $this_number
	cd $this_number

	# subset bed file to only include the gene we're interested in right now
	tail -n +$(( this_number + 1 )) ../Gene_Regions.bed | head -n 1 > ./Gene_Regions_no_header.bed

	# now subset the big vcf file so that it only includes genomic regions inside the gene of interest
	bcftools view --regions-file ./Gene_Regions_no_header.bed ../genotypes.MAF01.vcf.gz -O z -o ./genotypes.MAF01.GeneRegions.vcf.gz
	
	# decompress it
	gunzip ./genotypes.MAF01.GeneRegions.vcf.gz

	# process the vcf file with our awk script
	awk -f ../003_process_vcf.awk ./genotypes.MAF01.GeneRegions.vcf > ./genotypes.MAF01.GeneRegions.Parsed.txt

	# calculate p values for all possible qtls
	Rscript ../003_Find_sQTLs.R
	
	# move back up before creating a new folder
	cd ..
done
```

This will take a super long time to run. If you have access to a cluster and you can submit jobs using `qsub`, you can use the following version of the script to run each iteration of the loop as a separate job:

```bash
#!/bin/bash
#$ -N sQTL
#$ -t 1-4885

# save SGE_TASK_ID with a more user-friendly name
this_number=${SGE_TASK_ID}

# make a new directory and move to it
mkdir $this_number
cd $this_number

# subset bed file to only include the gene we're interested in right now
tail -n +$(( this_number + 1 )) ../Gene_Regions.bed | head -n 1 > ./Gene_Regions_no_header.bed

# now subset the big vcf file so that it only includes genomic regions inside the gene of interest
bcftools view --regions-file ./Gene_Regions_no_header.bed ../genotypes.MAF01.vcf.gz -O z -o ./genotypes.MAF01.GeneRegions.vcf.gz

# decompress it
gunzip ./genotypes.MAF01.GeneRegions.vcf.gz

# process the vcf file with our awk script
awk -f ../003_process_vcf.awk ./genotypes.MAF01.GeneRegions.vcf > ./genotypes.MAF01.GeneRegions.Parsed.txt

# calculate p values for all possible qtls
Rscript ../003_Find_sQTLs.R 

```

## Awk script

This script takes a VCF file and creates another file that tells us the number of copies of each potential QTL present in each subject. The output is a table similar to this:

|                  | Subject 1 | Subject 2 | Subject 3 | Subject 4 | &nbsp;&nbsp; ... &nbsp;&nbsp; |
|------------------|-----------|-----------|-----------|-----------|-------------------------------|
| Potential sQTL 1 | 0         | 1         | 0         | 1         | &nbsp;&nbsp; ... &nbsp;&nbsp; |
| Potential sQTL 2 | 1         | 0         | 0         | 2         | &nbsp;&nbsp; ... &nbsp;&nbsp; |
| Potential sQTL 3 | 0         | 2         | 1         | 0         | &nbsp;&nbsp; ... &nbsp;&nbsp; |

The number (between 0 and 2) indicates how many copies each subject has of a particular variant. This is the script I used:

```awk
BEGIN{OFS="\t";FS="\t"}
{
    l=$1"\t"$2"\t"$3;
    for(i=10;i<=NF;i++){
    if($1=="#CHROM"){
        g=$i
    } else {
        g="NA"
        if($i~"0/0:"){g=0}
        if($i~"1/0:"){g=1}
        if($i~"0/1:"){g=1}
        if($i~"1/1:"){g=2}
    }
    l=l"\t"g
    }
    print l;
}
```

## R script

This script looks at all the potential sQTL's for a splicing event of interest and performs a statistical test to determine whether they influence the behaviour of said splicing event. First, load the libraries we will need:


```r
library(lme4)
library(data.table)
```

Load the files generated in step 002 (`PSI.Estimates` and `GTEX.Dataset`):

```r
load("../Filtered.Annotations.PSI.Estimates.RData")
```
Load the BED file containing information about the region where this splicing event is located:

```r
Splicing.Events <- fread(input = "Gene_Regions_no_header.bed")
colnames(Splicing.Events) <- c("Chr", "Start", "End", "ID", "Gene")
```
Load the file produced by our awk script into an object called `Genotypes`:

```r
Genotypes <- fread("genotypes.MAF01.GeneRegions.Parsed.txt",
                   skip = "#CHROM")
```
Subset `Genotypes` so it only contains the subjects that passed the filter in step 002:

```r
Subjects.We.Can.Look.At <- which(colnames(Genotypes) %in% GTEX.Dataset$SUBJID)
Genotypes <- Genotypes[,.SD,,.SDcols=c(1,2,3,Subjects.We.Can.Look.At)]
```
Build a vector of the tissues we have samples for:

```r
Vector.of.Tissues.Analysed <- unique(as.character(GTEX.Dataset$SMTS))
```
Extract information relative to this particular exon splicing event:

```r
# some definitions for this particular event
This.Splicing.Event <- as.character(Splicing.Events$ID)
This.Event.Chromosome <- Splicing.Events$Chr
This.Event.Start <- Splicing.Events$Start
This.Event.End <- Splicing.Events$End
This.Event.PSI <- as.numeric(PSI.Estimates[This.Splicing.Event,])

# potential sQTLs associated with this splicing event
Potential.QTLs <- Genotypes[which(Genotypes$`#CHROM` == This.Event.Chromosome),]
Potential.QTLs <- Potential.QTLs[which(Potential.QTLs$POS < This.Event.End & Potential.QTLs$POS >= This.Event.Start), ]
```
For each variant in this gene, and for each tissue, I checked whether there were enough samples. Because I'm interested in checking the quantitative effect of the variant on exon inclusion, I want more than 10 samples without any copy of the variant, more than 10 samples with just one copy of the variant, and more than 10 samples with two copies of the variant.


PICTURE


If there were enough samples, built two generalised linear models with a binomial distribution and a random effect term (the subject ID):

<p align="center">
PSI ~ <b>Genotype</b> + Ischaemic Time + (1|Subject ID)
</p>

<p align="center">
PSI ~ Ischaemic Time + (1|Subject ID)
</p>


I then performed a likelihood ratio test between the two models to check whether adding the `Genotype` term significantly improves the model prediction.


```r
# stop if there are no potential sQTLs in this gene (super unlikely but just in case...)
if (nrow(Potential.QTLs) != 0) {
  # build a df for this particular exon splicing event
  This.Splicing.Event.Table <- GTEX.Dataset.Annotations
  This.Splicing.Event.Table$PSI <- This.Event.PSI
  
  # start an empty list
  sQTL.List <- vector(mode = "list", length = nrow(Potential.QTLs))
  
  # loop through all the potential QTLs and perform a statistical test if appropriate
  for (j in 1:nrow(Potential.QTLs)){
    This.Splicing.Event.Genotypes <- c()
    for(each.id in This.Splicing.Event.Table$SUBJID) {
      if (each.id %in% colnames(Potential.QTLs)) {
        This.Splicing.Event.Genotypes <- c(This.Splicing.Event.Genotypes,
                                           as.numeric(Potential.QTLs[j,
                                                                     each.id,
                                                                     with=FALSE]))
      } else {
        This.Splicing.Event.Genotypes <- c(This.Splicing.Event.Genotypes, NA)
      }
      
    }
    This.Splicing.Event.Table$Genotype <- This.Splicing.Event.Genotypes
    
    # start an empty list
    Tissues.List <- vector(mode = "list",
                           length = length(Vector.of.Tissues.Analysed))
    names(Tissues.List) <- Vector.of.Tissues.Analysed
    
    # go through each tissue and fill in corresponding element in Tissues.List
    for (each.tissue in Vector.of.Tissues.Analysed) {
      
      # a sub-dataframe only with rows corresponding to this tissue
      SubDF <- This.Splicing.Event.Table[which(as.character(This.Splicing.Event.Table$SMTS) == each.tissue),]
      SubDF <- SubDF[complete.cases(SubDF),]
      
      # only do a statistical test if we have 0, 1 and 2's; and if the no. of 2's is > 10
      if (length(table(SubDF$Genotype)) == 3) {
        if (table(SubDF$Genotype)[1] > 10  && table(SubDF$Genotype)[2] > 10  && table(SubDF$Genotype)[3] > 10 ){
          lm.psi <- SubDF$PSI
          lm.psi <- round(lm.psi*100)
          lm.psi <- cbind(lm.psi, 100-lm.psi)
          lm.genotype <- as.numeric(as.character(SubDF$Genotype))
          lm.ischaemic <- SubDF$SMTSISCH/1000
          lm.subject <- as.factor(SubDF$SUBJID)
          
          
          res<-try(suppressMessages(glmer(lm.psi ~ lm.genotype + lm.ischaemic + (1|lm.subject),
                                          family=binomial,
                                          nAGQ=25)),
                   silent=TRUE)
          res0<-try(suppressMessages(glmer(lm.psi ~ lm.ischaemic + (1|lm.subject),
                                           family=binomial,
                                           nAGQ=25)),
                    silent=TRUE)
          
          
          # do a likelihood-ratio test
          # the test statistic is asymptotically chi-sq distributed,
          # so I'll take the chi-sq statistic
          # between both models and calculate the p value from that
          P.Value <- anova(res, res0)$"Pr(>Chisq)"[2]
          Estimate <- summary(res)$coefficients[2,1]
          
          Tissues.List[[each.tissue]] <- c(P.Value, Estimate)
        } else {
          Tissues.List[[each.tissue]] <- NA
        }
        
      } else {
        Tissues.List[[each.tissue]] <- NA
      }
      
    }
    
    sQTL.List[[j]] <- Tissues.List
  }
  
  
  # save the sQTL list
  save(sQTL.List, file = "sQTL.RData")
}
```