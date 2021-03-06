# Raw files

In this document, I explain the code found in [001a\_download\_raw\_files.sh](001a_download_raw_files.sh) and [001b\_process\_raw\_files.R](001b_process_raw_files.R). The code in this markdown file is written in bash or in R, as indicated in each section. Please note that the files produced by this document are very big (> 700Mb) and so have not been included in the repository. These instructions should be enough to reproduce the missing files, but please feel free to contact me if you are having trouble with this.

## 1. Choice of tissues

To compare exon inclusion in two different tissues, we chose the two tissues (Brain and Skin) with the highest number of samples in the [GTEx Project](https://www.gtexportal.org/home/).

## 2. Download files

I downloaded 2 files from the GTEx Project [downloads page](https://www.gtexportal.org/gtex_analysis_v7/datasets):

* Sample annotations file: `GTEx_v7_Annotations_SampleAttributesDS.txt`
* Junction read counts file: `GTEx_Analysis_2016-01-15_v7_STARv2.4.2a_junctions.gct.gz`

To download these files, you'll need to create a GTEx login account for free.

## 3. Process sample annotations file

I used the sample annotations file to build a table (saved in a file called `IDs_Tissues.txt`) linking each GTEx sample ID with the tissue it comes from:

```bash
# bash script
cut -f1,6 GTEx_v7_Annotations_SampleAttributesDS.txt | sed s/-/\./g > IDs_Tissues.txt
```


## 4. Process junction read counts file

Before anything else, the junction read counts file has to be decompressed:

```bash
# bash script
gunzip GTEx_Analysis_2016-01-15_v7_STARv2.4.2a_junctions.gct.gz
```

To process the GTEx junction read counts file and extract PSI information about alternative splicing events throughout the genome, I used the Psichomics library in R (available on [GitHub](https://github.com/nuno-agostinho/psichomics) and [Bioconductor](https://bioconductor.org/packages/release/bioc/html/psichomics.html)). **What follows in the remainder of this document is R code.**

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
As well as the type of splicing event to quantify. Psichomics can quantify different types of [splicing events](http://rstudio-pubs-static.s3.amazonaws.com/359093_7f4afce0df5d48ba99eb0c05a9af8a00.html#quantifying-alternative-splicing), but in this case I was only interested in alternative splice site usage. For alternative 5' splice sites, you need to use `getSplicingEventTypes()[3]` as shown in the code below. For 3' splice sites, you need to use `getSplicingEventTypes()[4]`:

```r
# Events to quantify (just look at alt 5' splice sites for now)
eventType <- getSplicingEventTypes()[3]
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
PSU.Estimates <- quantifySplicing(annotation = annotation,
                                  junctionQuant = junctionQuant,
                                  eventType = eventType, 
                                  minReads = minReads)
```
**Warning**: because the junction read counts file is so big (~9Gb), the splicing quantification step can take a long time and your R session might run out of RAM. If that occurs, the easiest workaround is to split the junction file into 5-10 different files and process each of them separately, and merge the final processed files later.

## 5. Subset annotations data frame and save

GTEx did not count junction reads for all its sample. Therefore, I subsetted the annotations data frame so it only contained samples found in the `PSU.Estimates` data frame. To do this, I first loaded `IDs_Tissues.txt` into R:

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
# the PSU estimates dataframe
Common.Samples <- intersect(as.character(rownames(GTEX.Dataset)),
                            as.character(colnames(PSU.Estimates)))

# subset
GTEX.Dataset <- GTEX.Dataset[Common.Samples,]
PSU.Estimates <- PSU.Estimates[,Common.Samples]

```
And finally, save both R objects:

```r
save(GTEX.Dataset, PSU.Estimates, file = "Tissues_Compared_Datasets.RData")
```
