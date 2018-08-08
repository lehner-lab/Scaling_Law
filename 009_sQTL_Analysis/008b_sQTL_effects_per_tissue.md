# R code for sQTL effect per tissue


This document explains the code found inside [008b\_sQTL\_effects\_per\_tissue.R](008b_sQTL_effects_per_tissue.R). This script calculates the effect of a given sQTL in different tissues. Unless stated otherwise, all the code in this document is written in R.


## 1. Preparation

Before extracting sQTL effects, we need to prepare a table
Table with information about this sQTL + this splicing event including sample tissue, sample genotype and exon PSI. To do this, we first load the data table library:

```r
# load data handling lib
library(data.table)
```
If you remember the [previous document](008a_all_sQTL_effects_per_tissue.md), we called [008b\_sQTL\_effects\_per\_tissue.R](008b_sQTL_effects_per_tissue.R) with the following **bash** command: `Rscript ../../008b_sQTL_effects_per_tissue.R $sQTL_id`. The `$sQTL_id` argument contains the ID of the sQTL whose effect we want to check. We will pass this on to R using the `commandArgs` function:

```r
# sQTL ID is fed from the command line
args = commandArgs(trailingOnly=TRUE)
sQTL.ID <- args[1]
```
We have the sQTL ID, now we extract the splicing event ID:

```r
# splicing event ID taken from the BED file
Splicing.Event.ID <- as.character(read.table("Gene_Regions.bed")[1,4])
```
We'll now load `Annotations_PSIestimates.RData` in order to obtain information about the exon's PSI in each sample:

```r
# load GTEX.Dataset.Annotations and PSI.Estimates
load("../../Data/Annotations_PSIestimates.RData")
```
Before getting the PSI information, we need the index of this splicing event:

```r
# import a vector with the IDs of all 35153 splicing events analysed
All.Splicing.Events <- as.character(read.table("../../Gene_Regions.bed", header = T)$ID)

# which row contains our splicing event of interest?
Splicing.Event.Row <- which(All.Splicing.Events == Splicing.Event.ID)
```
Use the index obtained to extract all PSI values associated with our splicing events from `PSI.Estimates` and save them inside a vector called `This.Event.PSI`:

```r
# PSI values for this splicing event across all our samples
This.Event.PSI <- as.numeric(PSI.Estimates[Splicing.Event.Row,])
```
Now, we'll build a table for this splicing event containing:

* stuff from `GTEX.Dataset.Annotations`
* PSI
* Subject ID

```r
# Build a table with information about this splicing event
This.Splicing.Event.Table <- GTEX.Dataset.Annotations
This.Splicing.Event.Table$PSI <- This.Event.PSI
This.Splicing.Event.Table$Subject.ID <- gsub(pattern = "\\.",
                                             replacement = "-",
                                             x = This.Splicing.Event.Table$Subject.ID)
```
We will add a new column to this table containing information about the sQTL genotype in each sample. To do this, first load the parsed genotypes file:

```r
# load parsed genotypes file
Genotypes <- fread(input = "genotypes.MAF01.GeneRegions.Parsed.txt",
                   skip = "#CHROM")
```
Next, take all the sQTL genotype data for all the samples:
                 
```r
# row for this sQTL
sQTL.Row <- which(as.character(Genotypes$ID) == sQTL.ID)

# extract genotypes for this sQTL
This.Splicing.Event.Genotypes <- as.numeric(Genotypes[sQTL.Row,4:ncol(Genotypes)])
names(This.Splicing.Event.Genotypes) <- colnames(Genotypes)[4:ncol(Genotypes)]
```
Before adding it to `This.Splicing.Event.Table`, make sure that the sQTL vector is in the same order as the samples in the table itself:

```r
# empty vector which we'll fill in with genotype data, in
# the same order as the samples from This.Splicing.Event.Table
Genotypes.For.Table <- c()

# fill in the empty vector with sQTL genotype in the correct order
for(each.id in This.Splicing.Event.Table$Subject.ID) {
  if (each.id %in% names(This.Splicing.Event.Genotypes)) {
    Genotypes.For.Table <- c(Genotypes.For.Table, as.numeric(This.Splicing.Event.Genotypes[each.id]))
  } else {
    Genotypes.For.Table <- c(Genotypes.For.Table, NA)
  }
}
```
Add the sQTL genotype column to the table:

```r
This.Splicing.Event.Table$Genotype <- Genotypes.For.Table
```


## 2. Effect of sQTL in each tissue

Now that `This.Splicing.Event.Table` is complete, we can calculate the effect of all sQTLs in each tissue. First, build a vector with the names of all the different tissues we have data for:

```r
# get all possible tissues
All.Tissues <- unique(GTEX.Dataset.Annotations$Tissue)
```

Now, I'll initialise two empty vectors:

1. `Vector.Tissue.PSI.At.Gt.0` , which I'll fill in later with information about the exon's PSI in each tissue in the absence of the sQTL.
2. `Vector.Tissue.Effect.Of.QTL`, which I'll fill in with the sQTL effect in each tissue. 

```r
# initialise empty vector for PSI at genotype 0
Vector.Tissue.PSI.At.Gt.0 <- rep(NA, length(All.Tissues))
names(Vector.Tissue.PSI.At.Gt.0) <- as.character(All.Tissues)

# empty vector for sQTL effect
Vector.Tissue.Effect.Of.QTL <- rep(NA, length(All.Tissues))
names(Vector.Tissue.Effect.Of.QTL) <- as.character(All.Tissues)
```
We'll now iterate through the different tissues and fill in the vectors initialised above:

```r
# loop through the different tissues and get:
# - PSI of splicing event at gt == 0
# - Effect of sQTL per unit gt in all tissues
for (each.tissue in All.Tissues) {
  SubDF <- This.Splicing.Event.Table[which(This.Splicing.Event.Table$Tissue == each.tissue),]
  
  Which.Genotypes.To.Consider <- which(table(SubDF$Genotype) >= 7)
  
  if (length(Which.Genotypes.To.Consider) >= 2) {
    Genotype.To.Ignore <- which(table(SubDF$Genotype) < 7) - 1
    
    if (! length(Genotype.To.Ignore) == 0){
      SubDF <- SubDF[- which(SubDF$Genotype == Genotype.To.Ignore),]
    }
    
    PSI.At.Gt.0 <- lm(PSI ~ Genotype, data = SubDF)$coefficients[1]
    Vector.Tissue.PSI.At.Gt.0[each.tissue] <- PSI.At.Gt.0
    
    Effect.Of.QTL <- lm(PSI ~ Genotype, data = SubDF)$coefficients[2]
    Vector.Tissue.Effect.Of.QTL[each.tissue] <- Effect.Of.QTL
  }
  
}
```
Combine the two vectors, `Vector.Tissue.PSI.At.Gt.0` and `Vector.Tissue.Effect.Of.QTL ` into one long vector:

```r
# join everything into one long vector
This.sQTL.Information <- c(Vector.Tissue.PSI.At.Gt.0,
                           Vector.Tissue.Effect.Of.QTL)
```
Finally, save this vector:

```r
# save the vector
save(This.sQTL.Information, file = paste(paste("Effect.Of.sQTL.In.Different.Tissues",
                                               sQTL.ID,
                                               sep = "@"),
                                         ".RData",
                                         sep = ""))
```