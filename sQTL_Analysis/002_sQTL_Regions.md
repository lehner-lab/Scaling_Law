
# Defining regions where to hunt for sQTLs

In theory, a DNA variant could potentially be an sQTL for an exon anywhere else in the genome (even in a different chromosome). However, to limit the computational burden of the analysis, I only considered variants within the same gene as the exon of interest.

In this document I:

* identify suitable exons for the sQTL analysis
* determine the genomic regions I want to focus my analysis on

## Download genome annotations

In addition to the files downloaded in the previous document, I will also use the GENCODE annotations for the human genome release 19 (GRCh37.p13) (`gencode.v19.annotation.gff3.gz`, which can be downloaded from [this page](https://www.gencodegenes.org/releases/19.html)).

## Variable exons

To see whether the behaviour of sQTL's also scales depending on the starting PSI of the exon they affect, we need to find exons with very different levels of inclusion in different tissues.

I first built a table linking each GTEx sample to the tissue it comes from (as well as the ischaemic time of the sample, but this is for a downstream analysis later on). The following bash command does this:

```bash
cut -f1,6,9 GTEx_v7_Annotations_SampleAttributesDS.txt | sed s/-/\./g > IDs_Tissues_IschaemicTime.txt
```
I loaded `PSI.Estimates.RData` into R:

```r
load("PSI.Estimates.RData")
```

As well as the ID &rarr; Tissue table we just built:

```r
# build a dictionary to match GTEx IDs to tissue of origin
Dictionary <- read.table(file = "IDs_Tissues_IschaemicTime.txt",
                         header = TRUE,
                         sep = "\t")
                         
# only need the first 2 columns
Dictionary <- Dictionary[,1:2]

# rename columns, rows and convert from factor to character
colnames(Dictionary) <- c("ID", "Tissue")
Dictionary$ID <- as.character(Dictionary$ID)
Dictionary$Tissue <- as.character(Dictionary$Tissue)
rownames(Dictionary) <- Dictionary$ID
```
Next, I split the `PSI.Estimates` data frame into many data frames according to the sample tissue:

Vector with sample IDs from PSI object, and another with corresponding vectors

```r
# build a vector with the IDs of all samples in PSI.Estimates
IDs <- as.character(colnames(PSI.Estimates))
IDs <- IDs[2:length(IDs)]

# build a vector with the tissue of each sample in PSI.Estimates
Tissues <- Dictionary[IDs, "Tissue"]

# create an empty list with one element for each tissue
Split.PSI.Tables <- vector(mode = "list",
                           length = length(unique(Tissues)))
names(Split.PSI.Tables) <- unique(Tissues)

# fill the empty list with Build a PSI.Estimates table for each tissue
for (Each.Tissue in unique(Tissues)) {
  Columns <- which(Tissues == Each.Tissue)
  Split.PSI.Tables[[Each.Tissue]] <- PSI.Estimates[,c(1,(Columns+1))]
}

```
For each exon skipping event, calculate the median PSI in each tissue:

```r
# estimate median PSI per tissue
Median.PSIs.Per.Tissue <- lapply(X = Split.PSI.Tables,
                                 FUN = function(x){
                                   apply(X = x[,2:ncol(x)],
                                         MARGIN = 1,
                                         FUN = median,
                                         na.rm = TRUE)
                                   })

# matrix with all splicing events (rows) against tissue (columns)
Medians.Matrix <- do.call(cbind, Median.PSIs.Per.Tissue)
rownames(Medians.Matrix) <- as.character(rownames(PSI.Estimates))
```

As mentioned above, I am interested in exon splicing events with large PSI differences between tissues. To do this, I'll hunt for events where the PSI of the exon is very different from the rest in at least one tissue. However, I don't want this tissue to be associated with a specific sex (the differences in inclusion could then have to do with the sex of the subject). Furthermore, I also don't want the tissue to have fewer than 100 samples (i.e. I want to be really sure that the differences in tissues are real and not just due to a small sample size).

Therefore, before I proceed, I'll remove some tissues from `Medians.Matrix`:

```r
# define tissues I want to exclude
Bad.Tissues <- c("Ovary",
                 "Uterus",
                 "Vagina",
                 "Breast",
                 "Salivary Gland",
                 "Testis",
                 "Kidney",
                 "Fallopian Tube",
                 "Bladder",
                 "Cervix Uteri",
                 "Prostate")

# remove tissues I want to exclude
Medians.Matrix <- Medians.Matrix[,which(! colnames(Medians.Matrix) %in% Bad.Tissues)]
```

Now I can see which exon skipping events are the most variable. I defined 'variable exons' as those whose PSI was intermediate (between 40% and 60%) in at least one tissue, and either high (above 80%) or low (under 20%) in at least another tissue.

```r
# function to decide classify a splicing event according to the PSI
Categorise <- function(psi){
  Result <- NA
  if (! is.na(psi)){
    if (psi < 20){
      Result <- "Low"
    } else if (psi < 40){
      Result <- "MidLow"
    } else if (psi < 60){
      Result <- "Mid"
    } else if (psi < 80){
      Result <- "MidHigh"
    } else if (psi <= 100){
      Result <- "High"
    }
  }
  Result
}

# convert Medians.Matrix (with PSI numbers) into Categories.Matrix (with
# PSI categories instead)
Categories.Matrix <- apply(X = Medians.Matrix*100,
                           MARGIN = 1:2,
                           FUN = Categorise)
colnames(Categories.Matrix) <- colnames(Medians.Matrix)
rownames(Categories.Matrix) <- rownames(Medians.Matrix)

# use Categories.Matrix to see whether a given row (a splicing event)
# displays enough variability
Enough.Variability <- apply(X = Categories.Matrix,
                            MARGIN = 1,
                            FUN = function(x){
                              Result <- FALSE
                              if ("High" %in% x & "Mid" %in% x) {
                                Result <- TRUE
                              } else if ("Low" %in% x & "Mid" %in% x){
                                Result <- TRUE
                              }
                              Result
                            })

# so, which exons are variable?
Variable.Exons <- which(Enough.Variability)
Variable.Exons <- as.character(rownames(PSI.Estimates))[Variable.Exons]
```
Save `Variable.Exons` for later use. This vector contains the Psichomics IDs of all the splicing events I wanted to look at.

```r
# save
save(Variable.Exons, file = "Variable.Exons.RData")
```


## Exonic regions

Originally, I thought of only looking at variants located in the same exon whose PSI they affect. To do this, I loaded `Variable.Exons.RData` into R:

```r
load("Variable.Exons.RData")
```
I then used the ID of each splicing event to extract the chromosome where the event takes place, the start and end coordinates of the exon, as well as the name of the gene where the event occurs:

```r
# which chromosome
Chromosomes <- sapply(as.character(Variable.Exons),
                      function(x){
                        strsplit(x,"_")[[1]][2]
                      })

# starting coordinates
Start.BED <- sapply(as.character(Variable.Exons),
                    function(x){
                      strand <- strsplit(x,"_")[[1]][3]
                      if (strand == "+") {
                        start <- as.numeric(strsplit(x,"_")[[1]][5])
                      } else {
                        start <- as.numeric(strsplit(x,"_")[[1]][6])
                      }
                      start
                    })

# end coordinates
End.BED <- sapply(as.character(Variable.Exons),
                  function(x){
                    strand <- strsplit(x,"_")[[1]][3]
                    if (strand == "+") {
                      end <- as.numeric(strsplit(x,"_")[[1]][6])
                    } else {
                      end <- as.numeric(strsplit(x,"_")[[1]][5])
                    }
                    end
                  })

# gene symbols
Genes <- sapply(as.character(Variable.Exons),
                function(x){
                  split_id <- strsplit(x,"_")[[1]]
                  split_id[length(split_id)]
                })
```
And finally, save this all as a BED file for later use:

```r
# join everything in one data frame
Splicing.Events.BED <- data.frame(Chr = Chromosomes,
                                  Start = Start.BED,
                                  End = End.BED,
                                  ID = as.character(Variable.Exons),
                                  Gene = Genes)

# save as a BED file
write.table(x = Splicing.Events.BED,
            file = "Exonic_Regions.bed",
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")
```

## Gene regions

Not many variants (from the VCF file) were located in the regions that I just defined (i.e. inside the same exon they potentially regulate). I decided to amplify the region where I looked at potential sQTL's to include any variant located within the same gene as the exon splicing event. For this, I needed to use the human genome annotations I downloaded previously:

```bash
gunzip gencode.v19.annotation.gff3.gz
```
Take the coordinates corresponding to the start and end of the gene each of our exon splicing events are found in.

```bash
# an array with all the genes we are interested in
GENES=$(tail -n +2 Exonic_Regions.bed | cut -f 5)

c=0
for each_gene in $GENES; do
	# find the line corresponding to this gene
	each_gene_limits=$(grep "gene_name=${each_gene};" gencode.v19.annotation.gff3 | head -n 1 | cut -f 4,5)
	
	# take the start and the end coordinates for the whole gene
	gene_start=$(echo $each_gene_limits | awk '{print $1 - 1}')
	gene_end=$(echo $each_gene_limits | awk '{print $2}')
	
	# save these coordinates in an array
	all_starts[$c]=$gene_start
	all_ends[$c]=$gene_end
	
	# count +1
	c=$(($c+1))
done
```
Finally, create a new BED file with the new, updated coordinates:

```bash
# create the new bed file with gene regions
head -n 1 Exonic_Regions.bed > Gene_Regions.bed
paste <(tail -n +2 Exonic_Regions.bed) <(printf "%s\n" "${all_starts[@]}") <(printf "%s\n" "${all_ends[@]}") | awk 'BEGIN {FS="\t"} {print $1 FS $6 FS $7 FS $4 FS $5}' >> Gene_Regions.bed
```

## Prepare datasets for sQTL tests (next document)

First, I read `IDs_Tissues_IschaemicTime.txt` into R and tweaked it a little bit to contain all the information I will need downstream:

```r
library(data.table)

# load table with sample information
GTEX.Dataset <- fread(input = "IDs_Tissues_IschaemicTime.txt")

# create a new column with the subject ID
GTEX.Dataset $SUBJID <- sapply(as.character(All.Samples$SAMPID),
                               function(x){
                                 paste(strsplit(x,"\\.")[[1]][1:2],
                                       sep = "-",
                                       collapse = "-")
                               })

# factorise subject ID and tissue
GTEX.Dataset$SUBJID <- as.character(GTEX.Dataset$SUBJID) # subject ID
GTEX.Dataset$SMTS <- as.factor(GTEX.Dataset$SMTS) # tissue

# convert ischaemic time to numeric
GTEX.Dataset$SMTSISCH <- sapply(GTEX.Dataset$SMTSISCH,
                                function(this.string){
                                  gsub(pattern = "^\\.",
                                       replacement = "-",
                                       x = this.string,
                                       perl = T)
                                })
GTEX.Dataset$SMTSISCH <- as.numeric(GTEX.Dataset$SMTSISCH)
```
Next, I loaded `PSI.Estimates`. This is a very big file containing information from many exon splicing events I am not interested in. Therefore, I used `Gene_Regions.bed` to reduce `PSI.Estimates` such that it only includes exons I wanted to analyse:

```r
# load data
Splicing.Events <- fread(input = "Gene_Regions.bed")
load("PSI.Estimates.RData")

# subset the massive dataframe to only have events we're interested in
PSI.Estimates <- PSI.Estimates[as.character(Splicing.Events$ID),]

# also subset 
IDs.To.Keep <- intersect(as.character(colnames(Merged.PSI)), as.character(GTEX.Dataset$SAMPID))
PSI.Estimates <- PSI.Estimates[,IDs.To.Keep]
```
Finally, subset `GTEX.Dataset` so it includes the same samples as `PSI.Estimates`:

```r
rownames(GTEX.Dataset) <- as.character(GTEX.Dataset$SAMPID)
GTEX.Dataset <- GTEX.Dataset[SAMPID %in% IDs.To.Keep]
```
And save objects:

```r
save(GTEX.Dataset,
     PSI.Estimates,
     file="Filtered.Annotations.PSI.Estimates.RData")

```
