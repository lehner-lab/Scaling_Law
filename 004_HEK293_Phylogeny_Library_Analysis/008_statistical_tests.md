# Statistical tests

In this document I explain the code found in [008\_statistical\_tests.R](./008_statistical_tests.R), where I perform some simple statistics on the dataset.

Unless stated otherwise, all the code in this document is written in R.


## 1. Genotypes different from the ancestor

First load the dataset we will need for this:

```r
# load Whole.Dataset
load("002_whole_dataset.RData")

# low variance subset of the data
Low.Variance <- Whole.Dataset[which(Whole.Dataset$SD < 10),]
```
For all low-variance genotypes, perform a Wilcoxon rank sum test to compare the genotype's 9 PSI estimates (from the 9 biological replicates) against the ancestral PSI (96%):

```r
# wilcoxon test for difference with 96% PSI
P.Different.From.Ancestor <- apply(X = Low.Variance[,2:10],
                                   MARGIN = 1,
                                   FUN = function(x){
                                     wilcox.test(x = x, mu = 96)$p.value
                                   })
```
FDR-correct the p values:

```r
# fdr-correct the p values
FDR.Different.From.Ancestor <- p.adjust(p = P.Different.From.Ancestor,
                                        method = "fdr")
```
And count how many genotypes are different from the ancestor:

```r
# count how many are different from ancestor
Number.Significant.Differences <- sum(FDR.Different.From.Ancestor < 0.05)
Total.Genotypes <- length(FDR.Different.From.Ancestor)

Number.Significant.Differences/Total.Genotypes * 100
```


## 2. Genotypes different from the human sequence

For all low-variance genotypes, perform a Wilcoxon rank sum test to compare the genotype's 9 PSI estimates (from the 9 biological replicates) against the human PSI (59.5%):

```r
# wilcoxon test for difference with 59.5% PSI
P.Different.From.Human <- apply(X = Low.Variance[,2:10],
                                MARGIN = 1,
                                FUN = function(x){
                                  wilcox.test(x = x, mu = 59.5)$p.value
                                })
```
FDR-correct the p values:

```r
# fdr-correct the p values
FDR.Different.From.Human <- p.adjust(p = P.Different.From.Human,
                                     method = "fdr")
```
And count how many genotypes are different from the human PSI:

```r
# count how many are different from human
Number.Significant.Differences <- sum(FDR.Different.From.Human < 0.05)
Total.Genotypes <- length(FDR.Different.From.Human)

Number.Significant.Differences/Total.Genotypes * 100
```


## 3. Behaviour of T19G

To test how often T19G promotes inclusion or skipping, we first load the Final PSI - Starting PSI data set into R:

```r
# load final PSI - starting PSI dataset
load("004_final_vs_starting_psi_list.RData")
```
We extract the T-19-G table from this data set and call it `T19G.Table`:

```r
# extract table for mutation T-19-G
T19G.Table <- Final.Vs.Starting.PSI.List[["T-19-G"]]

# factor to character
T19G.Table$Genotype.Final <- as.character(T19G.Table$Genotype.Final)
T19G.Table$Genotype.Starting <- as.character(T19G.Table$Genotype.Starting)

# change the ID of the ancestral genotype
T19G.Table$Genotype.Starting[which(T19G.Table$Genotype.Starting == "")] <- "ANCESTOR"
```
Rename the row names of `Whole.Dataset` to the genotype IDs:

```r
# set Whole.Dataset row names to be the genotype ID
rownames(Whole.Dataset) <- as.character(Whole.Dataset$Mutation.IDs)
rownames(Whole.Dataset)[which(Whole.Dataset$Mutation.IDs == "")] <- "ANCESTOR"
```
Use these new IDs to take the nine starting PSI replicates and the nine final PSI replicates:

```r
# take all 9 replicates of each starting and final PSI in T19G.Table
Starting.PSI <- Whole.Dataset[as.character(T19G.Table$Genotype.Starting),2:10]
Final.PSI <- Whole.Dataset[as.character(T19G.Table$Genotype.Final),2:10]
```
Make a Delta PSI data frame with information from the nine replicates:

```r
# make a data frame with the difference
Delta.PSI <- Final.PSI - Starting.PSI
```
Perform a Wilcoxon test on the Delta PSIs against a median of 0:

```r
# calculate P values
P.Values.Delta.PSI <- apply(X = Delta.PSI,
                            MARGIN = 1,
                            FUN = function(x){
                              wilcox.test(x = as.numeric(x), mu = 0)$p.value
                            })
```
FDR-adjust the p values obtained:

```r
# FDR-correct
FDR.Values.Delta.PSI <- p.adjust(p = P.Values.Delta.PSI, method = "fdr")
```
Count how often T19G promotes inclusion and how often it promotes skipping:

```r
# how many backgrounds promote inclusion/skipping?
sum(FDR.Values.Delta.PSI < 0.05 & T19G.Table$Delta.PSI > 0)
sum(FDR.Values.Delta.PSI < 0.05 & T19G.Table$Delta.PSI < 0)
```