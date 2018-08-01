# Combine splicing events with significant sQTLs

This document explains the code found in [007\_combine\_splicing\_events\_with\_significant\_sQTLs.R](007_combine_splicing_events_with_significant_sQTLs.R), where I generate a table of unique significant Splicing Event ID - sQTL ID combinations. All the code written in this document is written in R.

## 1. Preparation
First load the data table library:

```r
# libraries i need
library(data.table)
```
Generate a vector with the different tissues analysed in this data set:

```r
# Load data
load("Annotations_PSIestimates.RData")
Vector.of.Tissues.Analysed <- unique(as.character(GTEX.Dataset.Annotations$Tissue))
```

Create an empty data table with four columns:

1. Splicing event ID
2. sQTL candidate ID
3. Tissue where the test was performed
4. p value

```r
# empty data table to fill with test results from previous step
Final.Data.Table <- data.table(Splicing.Event.ID = character(),
                               sQTL.ID = character(),
                               Tissue = character(),
                               P.Value = numeric())
```
We're now going to loop through the 35153 splicing events for which we tried to find an sQTL and add the results of those tests to `Final.Data.Table`:

```r
for (i in 1:35153){
    # code here
}
```
Store in a string the location of three files which we will use now:

```r
# location of BED file for this splicing event
# contains chromosome, start, end, id, gene
Gene.Regions.File <- paste("sQTL_tests",
                           as.character(i),
                           "/Gene_Regions.bed",
                           sep = "",
                           collapse = "")

# location of parsed genotype file for this splicing event
Genotypes.File <- paste("sQTL_tests",
                        as.character(i),
                        "/genotypes.MAF01.GeneRegions.Parsed.txt",
                        sep = "",
                        collapse = "")

# location of sQTL p values file
File.Location <- paste("sQTL_tests",
                       as.character(i),
                       "/sQTL.RData",
                       sep = "",
                       collapse = "")
```
Load the BED file for this specific splicing event into R:

```r
# load BED file into R
Splicing.Events <- fread(Gene.Regions.File)
colnames(Splicing.Events) <- c("Chr", "Start", "End", "ID", "Gene")
```
And use it to extract the start and end coordinates of the gene associated with this splicing event:

```r
# Coordinates of gene associated to this splicing event
This.Event.Start <- Splicing.Events$Start
This.Event.End <- Splicing.Events$End
```
Next, load the parsed genotypes file:

```r
# load genotypes file into R
Genotypes <- fread(Genotypes.File,
                   skip = "#CHROM")

# make sure snp positions fall within the gene boundaries.
# because of how bed files are encoded, the last position should not be included
Genotypes <- Genotypes[which(Genotypes$POS < This.Event.End & Genotypes$POS >= This.Event.Start), ]
```
Now, before moving on and extracting the p values for each sQTL-exon-tissue combination, we'll create the other 3 vectors with which we will fill `Final.Data.Table`:

```r
# Final.Data.Table - Tissue column
Data.Frame.Tissue <- rep(Vector.of.Tissues.Analysed,
                         length(Genotypes$ID))

# Final.Data.Table - Splicing.Event.ID column
This.Event.ID <- as.character(Splicing.Events$ID)
Data.Frame.Splicing.ID <- rep(This.Event.ID,
                              length(Data.Frame.Tissue))

# Final.Data.Table - sQTL.ID column
Data.Frame.QTL.ID <- rep(as.character(Genotypes$ID),
                         each = length(Vector.of.Tissues.Analysed))
```


## 2. Gather information

Now that we're done with the preparatory work, we will collect all the p values calculated in []() and put them together in one large data table. Note that we'll only proceed if there is an `sQTL.RData` file for this splicing event. In most cases there will be, but if it wasn't produced (eg. if there weren't any potential sQTLs within the same gene as the exon of interest) we won't be considering this splicing event any further:

```r
# make sure the sQTL.RData exists; otherwise stop
if (file.exists(File.Location)){
    # code here
}
```
Following is the code that goes inside the `if` statement. First, load `sQTL.RData` into R:

```r
# load sQTL.RData containing sQTL.List
load(file = File.Location)

# Subset sQTL to only include IDs found in Genotypes.
# This step is necessary because I might have removed an ID
# earlier when discussing the issues with BED format etc
sQTL.List <- sQTL.List[Genotypes$ID]
```
Next, another `if` statement. This second statement checks that `sQTL.List` is not empty (although it shouldn't originally be empty, it might have become empty after filtering in the previous step:

```r
# make sure sQTL.List is not empty
if (length(sQTL.List) > 0) {
    # code here
}
```
The first thing we do once we've passed this second condition is to build a vector with the contents of the only column in `Final.Data.Table` that we haven't yet prepared: the one corresponding to the p values:

```r
# Final.Data.Table - P.Value column
Data.Frame.P.Value <- unlist(lapply(sQTL.List,
                                    function(x){
                                      unlist(lapply(x,
                                                    function(y){
                                                      if (is.na(y[1])){
                                                        result <- NA
                                                      } else {
                                                        result <- y[1]
                                                      }
                                                      result
                                                    }))
                                    }))
```
Now that we have the four columns ready, we create a data table with them:

```r
# build a data table with the four columns
Temporary.Data.Table <- data.table(Splicing.Event.ID = Data.Frame.Splicing.ID,
                                   sQTL.ID = Data.Frame.QTL.ID,
                                   Tissue = Data.Frame.Tissue,
                                   P.Value = Data.Frame.P.Value)
```
And append this data table to `Final.Data.Table`:

```r
Final.Data.Table <- rbindlist(l = list(Final.Data.Table, Temporary.Data.Table))
```



## 3. Unique combinations of significant sQTLs


All p values were Bonferroni corrected and a snp was considered an sQTL (i.e. its effect was considered significant) when the Bonferroni q value was < 0.05.

```r
# Bonferroni correction!!
Final.Data.Table <- Final.Data.Table[complete.cases(Final.Data.Table),]
Final.Data.Table$FDR.Value <- p.adjust(Final.Data.Table$P.Value, "bonferroni")
```
Make a table with only the significant rows:

```r
# subset of Final.Data.Table containing only significant rows
Significant.Final.Data.Table <- Final.Data.Table[which(Final.Data.Table$FDR.Value < 0.05),]
```
Now build a data frame with unique combinations of significant splicing event ID -  sQTL ID pairs:

```r
# unique Splicing Event ID - sQTL ID combinations
Unique.Combos <- unique(paste(Significant.Final.Data.Table$Splicing.Event.ID,
                              Significant.Final.Data.Table$sQTL.ID,
                              sep = "@"))


# store those combinations in a data frame
Split.Into.Splicing.QTL.IDs <- strsplit(Unique.Combos,
                                        split = "@")
Unique.Combos.DF <- data.frame(Splicing.Events = sapply(Split.Into.Splicing.QTL.IDs,
                                                        function(x) x[1]),
                               sQTL.ID = sapply(Split.Into.Splicing.QTL.IDs,
                                                function(x) x[2]))
```
Save `Unique.Combos.DF` in a txt file:

```r
write.table(x = Unique.Combos.DF,
            file = "Data/007_Significant_Splicing_Events_sQTLs_Combos.txt",
            col.names = FALSE,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)

```

`007_Significant_Splicing_Events_sQTLs_Combos.txt` corresponds to Supplementary Table 5 in the paper and has been uploaded to the `/Data` folder.