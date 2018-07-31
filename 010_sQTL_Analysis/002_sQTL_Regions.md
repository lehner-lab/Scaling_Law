
# Defining regions where to hunt for sQTLs

In theory, a DNA variant could potentially be an sQTL for an exon anywhere else in the genome (even in a different chromosome). However, to limit the computational burden of the analysis, I only considered variants within the same gene as the exon of interest.

In this document I determine the genomic regions I want to focus my analysis on and prepare the datasets for the sQTL testing that will come later.


## Download files

I downloaded 2 files from the GTEx Project [downloads page](https://www.gtexportal.org/gtex_analysis_v7/datasets). One is freely accessible to any user who has a GTEx login account:

* **_Sample annotations file_**: `GTEx_v7_Annotations_SampleAttributesDS.txt`

The other file is protected and you'll need to apply for access through [dbgap](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v7.p2):

* **_Subject annotations file_**: `phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt`

I will also use the GENCODE annotations for the human genome release 19 (GRCh37.p13), which can be downloaded from [this page](https://www.gencodegenes.org/releases/19.html)):

* **_Human genome annotations file_**: `gencode.v19.annotation.gff3.gz`

## Obtain list of genes

Since I wanted to restrict the sQTL search to the gene where the alternative exon is located, I first extracted the IDs of the relevant genes. I loaded `PSI.Estimates.RData` into R:

```r
load("PSI.Estimates.RData")
```
I then used the ID of each splicing event to extract the chromosome where the event takes place, as well as the name of the corresponding gene:

```r
# which chromosome
Chromosomes <- sapply(as.character(rownames(PSI.Estimates)),
                      function(x){
                        strsplit(x,"_")[[1]][2]
                      })

# gene symbols
Genes <- sapply(as.character(rownames(PSI.Estimates)),
                function(x){
                  split_id <- strsplit(x,"_")[[1]]
                  split_id[length(split_id)]
                })
```

I built a data frame containing the information I just extracted:

```r
# join everything in one data frame
Splicing.Events.Info <- data.frame(Chr = Chromosomes,
                                   ID = as.character(Variable.Exons),
                                   Gene = Genes)
```
And saved it as a table:

```r
# save as a BED file
write.table(x = Splicing.Events.Info,
            file = "Splicing.Events.Info.txt",
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")
```


## Gene regions

For each alternative exon, I used the GENCODE annotations to extract the coordinates of its gene and build a BED file to use further downstream. First, decompress the annotations file:

```bash
gunzip gencode.v19.annotation.gff3.gz
```
Next, remove commented lines from the file:

```bash
grep -v "^#" gencode.v19.annotation.gff3 > table.gencode.v19.annotation.gff3
```
After cleaning the annotations file, we're ready to build the BED file from within R. We first load the libraries necessary for this:

```r
library(data.table)
library(dplyr)
```
Load the annotations file into R:

```r
Gencode.Annotations <- fread(input = "table.gencode.v19.annotation.gff3")
colnames(Gencode.Annotations) <- c("Chromosome",
                                   "Source",
                                   "Type",
                                   "Start",
                                   "End",
                                   "Score",
                                   "Strand",
                                   "Phase",
                                   "Attributes")
```
I subsetted the data table so that only coordinates for gene elements labelled as "gene" were kept, and then extracted the gene symbol for each gene and added that as another column:

```r
Gencode.Annotations.Genes <- Gencode.Annotations[Type == "gene"]
Gencode.Annotations.Genes$Gene <- sapply(Gencode.Annotations.Genes$Attributes,
                                         function(x){
                                           strsplit(strsplit(x,";")[[1]][6], "=")[[1]][2]
                                         })
```
Next, load the `Splicing.Events.Info.txt` table that we saved in the previous section:

```r
Splicing.Events <- fread(input = "Splicing.Events.Info.txt")
```
Add two columns to `Splicing.Events`: one with the coordinates corresponding to the start of the gene, and one with the coordinates of the end of the gene:

```r
Splicing.Events$Start <- sapply(Splicing.Events$Gene,
                                function(x){
                                  row.number <- which(Gencode.Annotations.Genes$Gene == x)
                                  if (length(row.number) > 1) {
                                    row.number <- row.number[1]
                                  }
                                  if (! length(row.number) == 0){
                                    Start <- Gencode.Annotations.Genes$Start[row.number] - 1
                                  } else {
                                    Start <- NA
                                  }
                                  Start
                                })

Splicing.Events$End <- sapply(Splicing.Events$Gene,
                              function(x){
                                row.number <- which(Gencode.Annotations.Genes$Gene == x)
                                if (length(row.number) > 1) {
                                  row.number <- row.number[1]
                                }
                                if (! length(row.number) == 0){
                                  End <- Gencode.Annotations.Genes$End[row.number]
                                } else {
                                  End <- NA
                                }
                                End
                              })
```
Reorder the columns so that they follow the standard BED format:

```r
Splicing.Events <- Splicing.Events %>% select(Chr, Start, End, ID, Gene)
```
Remove rows (splicing events) for which I could not find gene coordinates. These events mostly correspond to splicing events or genes labelled as "hypothetical" by Psichomics.

```r
Splicing.Events <- Splicing.Events[complete.cases(Splicing.Events),]

```
And save as a BED file:

```r
write.table(x = Splicing.Events,
            file = "Gene_Regions.bed",
                   quote = FALSE,
                   row.names = FALSE,
                   col.names = TRUE,
                   sep = "\t")

```


## Prepare datasets for sQTL tests (next document)

I pre-processed the sample annotations file to keep the fields I was interested in (sample ID's, tissue, ischaemic time) and saved it as `IDs_Tissues_IschaemicTime.txt`:

```bash
cut -f1,6,9 GTEx_v7_Annotations_SampleAttributesDS.txt | awk 'BEGIN {FS="\t";OFS="\t"}; {gsub(/-/, ".", $1)}{print $0}' > IDs_Tissues_IschaemicTime.txt


```
Similarly, I processed the subject annotations file to keep the fields corresponding to the subject ID, sex and age, and saved it as `Subject_Sex_Age.txt`:

```bash
grep -v -e '^#' -e '^$' phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt | ghead -n -1 | cut -f 2,4,5 | sed s/-/\./ > Subject_Sex_Age.txt
```

I loaded both of these files into R:

```r
library(data.table)
All.Samples <- fread("IDs_Tissues_IschaemicTime.txt")
All.Subjects <- fread("Subject_Sex_Age.txt")
```
and combined them in a data table called `GTEX.Dataset.Annotations`:
 
```r
# create a column with subject id's
All.Samples$SUBJID <- sapply(as.character(All.Samples$SAMPID),
                             function(x){
                               paste(strsplit(x, "\\.")[[1]][1:2],
                                     sep = ".",
                                     collapse = ".")
                             })

# and merge both data tables using common column (SUBJID)
GTEX.Dataset.Annotations <- merge(All.Samples, All.Subjects)

# rename columns of merged data table
colnames(GTEX.Dataset.Annotations) <- c("Subject.ID",
                                        "Sample.ID",
                                        "Tissue",
                                        "Ischaemic.Time",
                                        "Sex",
                                        "Age")
```

Next, I loaded `PSI.Estimates`. For some of the splicing events in this R object, I could not find the corresponding genomic regions (see **Gene Regions** above). Therefore, I used `Gene_Regions.bed` to reduce `PSI.Estimates` such that it only includes exons I can analyse:

```r
# load data
Splicing.Events <- fread(input = "Gene_Regions.bed")
load("PSI.Estimates.RData")

# subset the massive dataframe to only have events we're interested in
PSI.Estimates <- PSI.Estimates[as.character(Splicing.Events$ID),]

# also subset 
IDs.To.Keep <- intersect(as.character(colnames(Merged.PSI)), as.character(GTEX.Dataset.Annotations$Sample.ID))
PSI.Estimates <- PSI.Estimates[,IDs.To.Keep]
```
Finally, subset `GTEX.Dataset` so it includes the same samples as `PSI.Estimates`:

```r
rownames(GTEX.Dataset.Annotations) <- as.character(GTEX.Dataset.Annotations$Sample.ID)
GTEX.Dataset.Annotations <- GTEX.Dataset.Annotations[Sample.ID %in% Sample.IDs.To.Keep]
```
And save objects:

```r
save(GTEX.Dataset.Annotations,
     PSI.Estimates,
     file="Filtered.Annotations.PSI.Estimates.RData")

```
