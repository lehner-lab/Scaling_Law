# GTEx annotations

In our sQTL analysis, a generalised linear model was built to predict exon PSI using the following variables as predictors:

* Genotype (the number of copies of the sQTL candidate)
* Sex (male/female)
* Age (age of subject)
* Ischaemic time (the time after death before the sample was processed)

Information about the genotype will be obtained from the file processed in the [first step](001_filter_genotypes_file.md) of the pipeline, but we still don't have any information about the other three variables.

In this document, I will describe the code from [005a_gtex_annotations.sh](005a_gtex_annotations.sh) and [005b\_gtex\_annotations.R](005b\_gtex\_annotations.R), where we obtain information about the Sex, Age and Ischaemic time of each sample. Scripts included in this document are written in either bash or R, as indicated in the corresponding sections.

## 1. GTEx files

**All code from this section is written in bash.**

Before doing anything else, we need to download and pre-process the relevant files from GTEx.

### 1.1. Download

I downloaded 2 files from the GTEx Project [downloads page](https://www.gtexportal.org/gtex_analysis_v7/datasets). One is freely accessible to any user who has a GTEx login account:

* **Sample annotations file**: `GTEx_v7_Annotations_SampleAttributesDS.txt`

The other file is protected and you'll need to apply for access through [dbgap](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v7.p2):

* **Subject annotations file**: `phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt`


### 1.2. Pre-process downloaded files

I pre-processed the sample annotations file to keep the fields I was interested in (sample ID's, tissue, ischaemic time) and saved it as `IDs_Tissues_IschaemicTime.txt`:

```bash
# pre-process sample annotations file
cut -f1,6,9 GTEx_v7_Annotations_SampleAttributesDS.txt | awk 'BEGIN {FS="\t";OFS="\t"}; {gsub(/-/, ".", $1)}{print $0}' > IDs_Tissues_IschaemicTime.txt
```
Similarly, I processed the subject annotations file to keep the fields corresponding to the subject ID, sex and age, and saved it as `Subject_Sex_Age.txt`:

```bash
# pre-process subject annotations file
grep -v -e '^#' -e '^$' phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt | head -n -1 | cut -f 2,4,5 | sed s/-/\./ > Subject_Sex_Age.txt
```

## 2. Build GTEx annotations R object

**All code in this section is written in R.**

Now that we have two txt files containing the information we want, we will combine them into one table. First, load the data table library in R:

```r
# data table package
library(data.table)
```
Read `IDs_Tissues_IschaemicTime.txt` and `Subject_Sex_Age.txt` into R:

```r
# load tables in R
All.Samples <- fread("IDs_Tissues_IschaemicTime.txt")
All.Subjects <- fread("Subject_Sex_Age.txt")
```
The subject ID consists of the first 2 components of the sample ID. In the `All.Samples` data table, create a column containing the subject ID:

```r
# subject ID column
All.Samples$SUBJID <- sapply(as.character(All.Samples$SAMPID),
                             function(x){
                               paste(strsplit(x, "\\.")[[1]][1:2],
                                     sep = ".",
                                     collapse = ".")
                             })
```
And now merge both tables into an object called `GTEX.Dataset.Annotations`:

```r
# note: for some reason no annotations for bone marrow subjects.
# Also ID has a different format. These samples are excluded in the merging step.
GTEX.Dataset.Annotations <- merge(All.Samples, All.Subjects)
colnames(GTEX.Dataset.Annotations) <- c("Subject.ID",
                                        "Sample.ID",
                                        "Tissue",
                                        "Ischaemic.Time",
                                        "Sex",
                                        "Age")

rownames(GTEX.Dataset.Annotations) <- as.character(GTEX.Dataset.Annotations$Sample.ID)
```


## 3. Filter GTEx annotations R object

**All code in this section is written in R.**

In the previous section, we built a table containing Sex, Age, Ischaemic Time and Tissue information for every sample in the GTEx data set. However, as explained in a [previous document](002_common_subjects.md), we are not working with the complete set of subjects because some were not genotyped and some do not have a junction read counts file. Therefore, we're now going to filter these "bad" samples out from `GTEX.Dataset.Annotations`:

```r
# load PSI.Estimates
load("PSI.Estimates.RData")

# What samples do we want to keep?
Sample.IDs.To.Keep <- intersect(as.character(colnames(PSI.Estimates)),
                                as.character(GTEX.Dataset.Annotations$Sample.ID))

# filter GTEX.Dataset.Annotations
GTEX.Dataset.Annotations <- GTEX.Dataset.Annotations[Sample.ID %in% Sample.IDs.To.Keep]
```
Finally, make sure rows and columns in `PSI.Estimates` match our other datasets (they should, but just in case):

```r
# make sure PSI.Estimates columns in the same order
PSI.Estimates <- PSI.Estimates[, as.character(GTEX.Dataset.Annotations$Sample.ID)]

# load splicing events table
Splicing.Events <- fread(input = "Gene_Regions.bed")

# make sure PSI.Estimates has the rows in the same order as Splicing.Events
PSI.Estimates <- PSI.Estimates[as.character(Splicing.Events$ID),]
```
And save everything:

```r
# save
save(GTEX.Dataset.Annotations,
     PSI.Estimates,
     file="Filtered_Annotations_PSIestimates.RData")

```