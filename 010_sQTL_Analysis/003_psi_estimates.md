# PSI estimates

In this document, I go over the code from [003\_psi\_estimates.R](003_psi_estimates.R), where I process the junction read counts file downloaded earlier to estimate the PSI of exon skipping events throughout the genome. All the code included in this document is written in R.

_____________________________________

To process the GTEx junction read counts file and extract PSI information about alternative splicing events throughout the genome, I used the Psichomics library in R (available on [GitHub](https://github.com/nuno-agostinho/psichomics) and [Bioconductor](https://bioconductor.org/packages/release/bioc/html/psichomics.html)).

First, load the library:

```r
library(psichomics)
```
Next, we have to specify what annotations we want to use:

```r
# Load Human (hg19/GRCh37 assembly) annotation
human <- listSplicingAnnotations()[[1]]
annotation <- loadAnnotation(human)
```
As well as the type of splicing event to quantify. Psichomics can quantify different types of [splicing events](http://rstudio-pubs-static.s3.amazonaws.com/359093_7f4afce0df5d48ba99eb0c05a9af8a00.html#quantifying-alternative-splicing), but in this case I was only interested in alternative exon usage:

```r
# Events to quantify (just look at alt exons for now)
eventType <- getSplicingEventTypes()[1]
```
Psichomics also requires the user to specify the minimum number of read counts needed for a splicing event in a given sample to be quantified. I set this to 10:

```r
# minimum read counts
minReads <- 10
```
Load the junction read counts file into R:

```r
# where is my file?
Junction.File.Path <- "GTEx_Analysis_2016-01-15_v7_STARv2.4.2a_junctions.gct"
  
# load it
GTEx.Junctions <- loadGtexData(junctionQuant = Junction.File.Path)
junctionQuant <- GTEx.Junctions$GTEx$"Junction quantification"
```
And quantify all the splicing events in the genome:

```r
# Quantify
PSI.Estimates <- quantifySplicing(annotation = annotation,
                                  junctionQuant = junctionQuant,
                                  eventType = eventType, 
                                  minReads = minReads)
```
**Warning**: because the junction read counts file is so big (~9Gb), the splicing quantification step can take a long time and your R session might run out of RAM. If that occurs, the easiest workaround is to split the junction file into 5-10 different files and process each of them separately, and merge the final processed files later.

Subset `PSI.Estimates` so that it only contains samples that are found both in the junctions and the genotypes files:

```r
# load vector with allowed subjects (common to genotypes + junctions file)
Allowed.Subjects <- as.character(read.table(file = "Common_Subjects.txt")$V1)

# all Samples from the PSI data frame
All.Samples <- colnames(PSI.Estimates)

# which of these samples can I work with?
Can.I.Work.With.This.Sample <- sapply(All.Samples,
                                      function(x){
                                        this.subject <- paste(strsplit(x = x,
                                                                       split = "\\.")[[1]][1:2],
                                                                       collapse = "-")
                                        this.subject %in% Allowed.Subjects
                                      })

Allowed.Samples <- which(Can.I.Work.With.This.Sample)

# update PSI.Estimates so it only contains samples we can work with
PSI.Estimates <- PSI.Estimates[,Allowed.Samples]
```
Some splicing events could not be quantified in any of the samples we are working with (i.e. the event returned an `NA` in every single sample). I removed these events from my dataset as they don't provide any useful information:

```r
# Look at mean NA value per splicing event
NAs.Per.Splicing.Event <- apply(X = PSI.Estimates,
                                MARGIN = 1,
                                FUN = function(x) mean(is.na(x)))

# select only events where I don't have any NAs
Not.All.NAs <- which(NAs.Per.Splicing.Event < 1)
PSI.Estimates <- PSI.Estimates[Not.All.NAs,]
```
And save the file.

```r
# save
save(PSI.Estimates, file = "PSI.Estimates.RData")
```
