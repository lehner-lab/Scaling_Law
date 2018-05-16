# Raw files


To compare exon inclusion in two different tissues, we chose the two tissues (Brain and Skin) with the highest number of samples in the [GTEx Project](https://www.gtexportal.org/home/). 

## Download files

I downloaded 2 files from the GTEx Project [downloads page](https://www.gtexportal.org/gtex_analysis_v7/datasets):

* Sample annotations file: `GTEx_v7_Annotations_SampleAttributesDS.txt`
* Junction read counts file: `GTEx_Analysis_2016-01-15_v7_STARv2.4.2a_junctions.gct.gz`

To download these files, you'll need to create a GTEx login account for free.

## Process sample annotations file

I used the sample annotations file to build a table linking each GTEx sample ID with the tissue it comes from:

```bash
# bash script
cut -f1,6 GTEx_v7_Annotations_SampleAttributesDS.txt | sed s/-/\./g > IDs_Tissues.txt
```


## Process junction read counts file

Before anything else, the junction read counts file has to be decompressed:

```bash
# bash script
gunzip GTEx_Analysis_2016-01-15_v7_STARv2.4.2a_junctions.gct.gz
```

To process the GTEx junction read counts file and extract PSI information about alternative splicing events throughout the genome, I used the Psichomics library in R (available on [GitHub](https://github.com/nuno-agostinho/psichomics) and [Bioconductor](https://bioconductor.org/packages/release/bioc/html/psichomics.html)). What follows in the remainder of this document is R code.

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

`PSI.Estimates` is an R data frame where columns represent different samples and rows represent different splicing events. I removed any splicing events (rows) which could not be quantified even once across all the samples:

```r
# Look at mean NA value per splicing event
NAs.Per.Splicing.Event <- apply(PSI.Estimates,
                                1,
                                function(x) mean(is.na(x)))

# select only events where I don't have any NAs
No.NAs <- which(NAs.Per.Splicing.Event == 0)
PSI.Estimates <- PSI.Estimates[No.NAs,]
```

## Subset annotations data frame and save

GTEx did not count junction reads for all its sample. Therefore, I subsetted the annotations data frame so it only contained samples found in the `PSI.Estimates` data frame. To do this, I first loaded `IDs_Tissues.txt` into R:

```r
# Load table with sample information
GTEX.Dataset <- read.table(file = "IDs_Tissues.txt",
                           header = T,
                           sep = "\t")
# rename columns & rows
colnames(GTEX.Dataset) <- c("ID", "Tissue")
rownames(GTEX.Dataset) <- as.character(GTEX.Dataset$ID)
```
And subset the annotations data frame:

```r
# find samples that are found in the annotations data frame and
# the PSI estimates dataframe
Common.Samples <- intersect(as.character(rownames(GTEX.Dataset)),
                            as.character(colnames(PSI.Estimates)))

# subset
GTEX.Dataset <- GTEX.Dataset[Common.Samples,]
PSI.Estimates <- PSI.Estimates[,Common.Samples]

```
And finally, save both R objects:

```r
save(GTEX.Annotations, PSI.Estimates, file = "Tissues_Compared_Datasets.RData")
```