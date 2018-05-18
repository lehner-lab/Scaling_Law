# Raw files


In this document, I define an sQTL as a QTL that alters the levels of inclusion of a particular exon. To analyse the dependence of sQTL effects on the starting level of exon inclusion, we chose to analyse data from the [GTEx Project](https://www.gtexportal.org/home/). This allowed us to look at the effect of the same sQTL in different tissues (i.e. different 'starting PSIs' if the exon studied is differentially included in different tissues).

Please note that many of the intermediate files produced in this document contain sensitive information extracted from a file whose access is restricted. Furthermore, many of these intermediate files are very big. Therefore, I have chosen to upload only the final processed files. The instructions here should be enough to reproduce the missing files, but please feel free to contact me if you are having trouble with this.

## Download files

I downloaded 3 files from the GTEx Project [downloads page](https://www.gtexportal.org/gtex_analysis_v7/datasets). Two of them are freely accessible to any user who has a GTEx login account:

* Sample annotations file: `GTEx_v7_Annotations_SampleAttributesDS.txt`
* Junction read counts file: `GTEx_Analysis_2016-01-15_v7_STARv2.4.2a_junctions.gct.gz`

The third file is protected and you'll need to apply for access through [dbgap](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v7.p2):

* Genotype matrix in VCF format: `phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1.GRU.tar` (the file we're interested in is `GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz` and is located inside this massive tar file)



## Obtain variants common in the population

The genotype calls file is huge. I subsetted to include only potential QTLs that are common in the population. The first thing I did was to rename the file and move it to the `Data` folder:

```bash
# rename and move to Data/ folder
mv GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz genotypes.vcf.gz
mv genotypes.vcf.gz Data/
```
Next, I subsetted the VCF file to include only potential QTLs that occur in at least 1% of all subjects. To do this, I used [BCFtools](https://samtools.github.io/bcftools/):

```bash
# subset the file to only include common snps, indels etc
bcftools view --min-af 0.01:minor Data/genotypes.vcf.gz -O z -o Data/genotypes.MAF01.vcf.gz  # 5-6 hours! from 120gb to 37gb
```
**Warning**: depending on how the VCF file was handled, you might get an error further downstream, saying something like *Could not load the index*. This error has been [reported before](https://github.com/samtools/bcftools/issues/129) and it is caused by your file having been compressed with gzip instead of [bgzip](http://www.htslib.org/download/). The fix is as follows:

```bash
# error reported in https://github.com/samtools/bcftools/issues/129
# the fix is this following line
zcat genotypes.MAF01.vcf.gz | bgzip -c > NEW.genotypes.MAF01.vcf.gz && tabix NEW.genotypes.MAF01.vcf.gz # file was compressed with gzip instead of bgzip ? this step took ~3-4 hours

# and rename
rm genotypes.MAF01.vcf.gz
mv NEW.genotypes.MAF01.vcf.gz genotypes.MAF01.vcf.gz
```


## Subjects in both the junctions and genotype files

Not all GTEx subjects had their genome sequenced and not all had their splicing junction reads counted. Therefore, to carry out an analysis involving both junctions and genotypes, I need to find out which subjects are common to both datasets.

```bash
# get subjects from junction read counts file
gunzip GTEx_Analysis_2016-01-15_v7_STARv2.4.2a_junctions.gct.gz
tail -n +3 GTEx_Analysis_2016-01-15_v7_STARv2.4.2a_junctions.gct | head -n 1 | awk 'BEGIN{OFS="\n"}{$1=$1; print $0}' | tail -n +3 | cut -d'-' -f 2 | awk '{print "GTEX-"$0}' | sort | uniq > Junctions_Subjects.txt

# get subjects from genotypes file
gunzip genotypes.MAF01.vcf.gz
grep -v "^##" genotypes.MAF01.vcf | head -n 1 | awk 'BEGIN{OFS="\n"}{$1=$1; print $0}' | tail -n +4 | sort | uniq > Genotypes_Subjects.txt

# intersect of both files
comm -12 Junctions_Subjects.txt Genotypes_Subjects.txt > Common_Subjects.txt
```


## Process junction read counts file

To process the GTEx junction read counts file and extract PSI information about alternative splicing events throughout the genome, I used the Psichomics library in R (available on [GitHub](https://github.com/nuno-agostinho/psichomics) and [Bioconductor](https://bioconductor.org/packages/release/bioc/html/psichomics.html)). What follows in the remainder of this section is R code.

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

