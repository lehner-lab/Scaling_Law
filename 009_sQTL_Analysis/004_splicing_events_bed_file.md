# Build a splicing events BED file

In this document I explain the code from [004a\_download\_gencode\_annotations.sh](004a_download_gencode_annotations.sh) and [004b\_splicing\_events\_bed\_file.R](004b_splicing_events_bed_file.R), where I build a table with the following format:

| Chr | Start     | End       | ID                                                           | Gene  |
|-----|-----------|-----------|--------------------------------------------------------------|-------|
| 10  | 101088855 | 101154087 | SE\_10\_+\_101124759\_101128375\_101128437\_101136812\_CNNM1 | CNNM1 |
| 10  | 101462314 | 101515891 | SE\_10\_+\_101507147\_101510126\_101510153\_101514286\_CUTC  | CUTC  |
| 10  | 102495359 | 102589698 | SE\_10\_+\_102506060\_102506772\_102507776\_102509503\_PAX2  | PAX2  |
| ... | ...       | ...       | ...                                                          | ...   |

Where ID refers to the splicing event ID and the start and end positions refer to the start and end positions of the gene where the splicing event happens.

Why do we need to do this? In theory, a DNA variant could potentially be an sQTL for an exon anywhere else in the genome (even in a different chromosome). However, to limit the computational burden of the analysis, I only considered variants within the same gene as the exon of interest.

The code described in this document is written in bash and in R, as stated in each section.



## 1. Download gene annotations

**The code in this section is written in bash.**

I used the GENCODE annotations for the human genome release 19 (GRCh37.p13), which can be downloaded from [this page](https://www.gencodegenes.org/releases/19.html)):

* **Human genome annotations file**: `gencode.v19.annotation.gff3.gz`

Download and extract the file:

```bash
# download file
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz

# unpack the file
gunzip gencode.v19.annotation.gff3.gz
```


## 2. Build a table with splicing events

**The code in this section is written in R.**

First, load the `PSI.Estimates` object created in an earlier step:

```r
# Load massive data frame with all the psi values
load("PSI.Estimates.RData")
```
For each splicing event, extract the chromosome it is found in:

```r
# chromosome id
Chromosomes <- sapply(as.character(rownames(PSI.Estimates)),
                      function(x){
                        strsplit(x,"_")[[1]][2]
                      })
```
as well as the gene to which it belongs:

```r
# gene symbol
Genes <- sapply(as.character(rownames(PSI.Estimates)),
                function(x){
                  split_id <- strsplit(x,"_")[[1]]
                  split_id[length(split_id)]
                })
```
Combine this information, along with the splicing event ID, in a data frame:

```r
# build a table with the following information:
# - chromosome (column 1)
# - id of splicing event (column 2)
# - gene symbol where splicing event is found
Splicing.Events.Table <- data.frame(Chr = Chromosomes,
                                    ID = as.character(rownames(PSI.Estimates)),
                                    Gene = Genes)
```


## 3. Add gene start/end positions to splicing events table

**The code in this section is written in R.**

To add the gene start and end positions to the data frame we have been building, we'll first load two data handling libraries:

```r
# some data handling libraries
library(data.table)
library(dplyr)
```
Load the GENCODE annotations file downloaded earlier into an object called `Gencode.Annotations`:

```r
# load Gencode annotations file using 'fread' function
Gencode.Annotations <- fread(input = "table.gencode.v19.annotation.gff3")

# add column names
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
Subset `Gencode.Annotations` to include only rows corresponding to gene features:

```r
# we're only interested in genes
Gencode.Annotations.Genes <- Gencode.Annotations[Type == "gene"]
```
Add a column to `Gencode.Annotations.Genes` containing the symbol of the gene in each row:

```r
# extract gene symbols
Gencode.Annotations.Genes$Gene <- sapply(Gencode.Annotations.Genes$Attributes,
                                         function(x){
                                           strsplit(strsplit(x,";")[[1]][6], "=")[[1]][2]
                                         })
```
Now we have two objects, `Gencode.Annotations.Genes` and `Splicing.Events.Table`, which both have a 'Gene' column containing the gene symbol/name. I'll use this common column to go through each gene symbol in `Splicing.Events.Table` and extract its start/end coordinates from `Gencode.Annotations.Genes`:

```r
# add gene start positions to Splicing.Events.Table
Splicing.Events.Table$Start <- sapply(Splicing.Events.Table$Gene,
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

# add gene end positions to Splicing.Events.Table
Splicing.Events.Table$End <- sapply(Splicing.Events.Table$Gene,
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
Rearrange the columns so that they fit the BED format:

```r
# shuffle Splicing.Events.Table columns around to fit BED format
Splicing.Events.Table <- Splicing.Events.Table %>% select(Chr, Start, End, ID, Gene)
Splicing.Events.Table <- Splicing.Events.Table[complete.cases(Splicing.Events.Table),]
```
And save the file:

```r
# save BED file
write.table(x = Splicing.Events.Table,
            file = "Gene_Regions.bed",
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")
```
