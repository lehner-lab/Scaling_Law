# Mathematical model of exon competition


## 1. Sigmoidal relationship between PSI and log k6

Equation 1 of our model relates exon 6's levels of inclusion with splicing efficiency parameters k6 and k6, as well as with a time delay parameter tau:

<p align="center">
  <img src="Equations/equation_1.gif">
</p>

The relationship between exon 6 PSI and fold-changes in k6 is sigmoidal. To show this, we'll plot PSI vs k6 while fixing k7 to 1 and tau to 0:

```r
# little function to represent equation 1 from the main text
Calculate.PSI <- function(k6, k7, tau){
  1 - exp(-1*k6*tau)*(k7/(k7+k6))
}

# fix tau and k7:
Fixed.Tau <- 0
Fixed.k7 <- 1

# a bunch of k6 values we'll test
Range.k6.To.Test <- 10^(seq(-3,3,0.005))

# the corresponding PSI values
Calculated.PSIs <- Calculate.PSI(k6 = Range.k6.To.Test,
                                 k7 = Fixed.k7,
                                 tau = Fixed.Tau)

# plot
plot(log2(Range.k6.To.Test),
     Calculated.PSIs*100,
     type = "l",
     col = "gray90",
     lwd = 3,
     xlab = "log2 k6",
     ylab = "PSI",
     las = 1)
```
<p align="center">
  <img width = 450 height = 450 src="Figures/001_sigmoid_psi_logK6.png">
  <br> Figure 4C
</p>



## 2. Final vs starting PSI

To visualise the relationship between final and starting PSIs, I used this equation:


<p align="center">
  <img src="Equations/Final_vs_Starting.gif">
</p>

to draw Final vs Starting PSI lines, each corresponding to a different value of A. To plot this, use the following code:

```r
# start an empty plot
par(pty = "s")
plot(NULL,
     xlim = c(0,100),
     ylim = c(0,100),
     xlab = "Starting PSI",
     ylab = "Final PSI",
     main = "Final vs Starting PSI",
     las = 1)

# draw lines
for (i in seq(-6,6,1)) {
  x.vector <- seq(0,1,0.001)
  y.vector <- (exp(i)*x.vector) / (1 + exp(i)*x.vector - x.vector)
  lines(100*x.vector,
        100*y.vector,
        col = "gray90",
        lwd = 3)
}
```

<p align="center">
  <img width = 450 height = 450 src="Figures/002_final_vs_starting_psi.png">
  <br> Figure 4D
</p>

## 3. Delta PSI vs starting PSI

This is very similar to what I did above. I now use the following equation instead:

<p align="center">
  <img src="Equations/Delta_vs_Starting.gif">
</p>

And the code to draw the graph is:

```r
# start an empty plot
par(pty = "s")
plot(NULL,
     xlim = c(0,100),
     ylim = c(-100,100),
     xlab = "Starting PSI",
     ylab = "Delta PSI",
     main = "Delta PSI and Starting PSI",
     las = 1)

# draw lines
for (i in seq(-6,6,1)) {
  x.vector <- seq(0,1,0.001)
  y.vector <- (exp(i)*x.vector) / (1 + exp(i)*x.vector - x.vector)
  y.vector <- y.vector - x.vector
  lines(100*x.vector,
        100*y.vector,
        col = "gray90",
        lwd = 3)}

```

<p align="center">
  <img width = 450 height = 450 src="Figures/003_delta_vs_starting_psi.png">
  <br> Figure 4D
</p>


## 4. Starting PSI at which maximum effect occurs

The relationship between the maximum effect size of a mutation and the starting PSI at which that effect is observed was derived to be:

<p align="center">
  <img src="Equations/Max_Effect.gif">
</p>

We can use this equation to visualise what this relationship looks like:

```r
# a bunch of maximum effect sizes
Max.Delta.PSIs <- seq(-100, 100, 20)

# the corresponding starting PSIs
Starting.PSIs <- 100 *(0.5 - 0.5 * Max.Delta.PSIs)

# plot
plot(NULL,
     ylim = c(0,100),
     xlim = c(-100,100),
     ylab = "Starting PSI",
     xlab = "Max Delta PSI")
par(new=T)
plot(y = Starting.PSIs,
     x = Max.Delta.PSIs,
     pch = 19,
     col = "gray80",
     axes = F,
     xlab = "",
     ylab = "",
     ylim = c(0,100),
     xlim = c(-100,100),
     cex = 1.5)
```

<p align="center">
  <img width = 450 height = 450 src="Figures/004_startingPSI_vs_max_effect.png">
  <br> Figure 4F
</p>


To see whether we observe this same relationship in our data set, we fitted non-parametric curves (principal curves) to the Final vs Starting PSI plots from Figure S3. We took the maximum effect size predicted by the non-parametric fit as an estimate of the mutation's maximum effect size and plotted the starting PSI where this effect is observed versus the effect itself. To do this, we first load the library needed for the fit, as well as the relevant datasets:

```r
# principal curves library
library(analogue)

# load datasets
load("../004_HEK293_Phylogeny_Library_Analysis/002_whole_dataset.RData")
load("../004_HEK293_Phylogeny_Library_Analysis/004_final_vs_starting_psi_list.RData")
```
A vector with the 12 mutations that occur in our library:

```r
# vector containing the identity of the 12 point mutations found in the library
Single.Mutations <- c("C-41-G",
                      "T-49-C",
                      "T-24-C",
                      "T-19-G",
                      "G-51-C",
                      "C-18-G",
                      "C-32-T",
                      "C-18-T",
                      "C-39-T",
                      "G-26-T",
                      "G-35-T",
                      "G-44-A")
```
Also make a note of which mutations have two clearly distinct behaviours. This is because we will fit a different curve to each of these behaviours:

```r
# mutations with two clearly distinct behaviours
Singles.With.Multiple.Paths <- c("C-18-G",
                                 "T-19-G",
                                 "C-32-T",
                                 "T-49-C",
                                 "G-51-C")

# epistatic partners of the above
Epistasis.Table <- data.frame(Main = Singles.With.Multiple.Paths,
                              Epistatic = c("T-19-G",
                                            "C-18-G",
                                            "G-35-T",
                                            "G-51-C",
                                            "T-49-C"))
```
Build a vector with the low-variance genotypes:

```r
# build a vector with low SD genotypes
Low.Noise.Dataset <- Whole.Dataset[which(Whole.Dataset$SD < 10),]
Low.SD.Variants <- as.character(Low.Noise.Dataset$Mutation.IDs)
```
Create an empty list called `Principal.Curves`, where we will store the non-parametric curves generated:

```r
Principal.Curves <- vector(mode = "list")
```
Now we are ready to generate the principal curves for all our mutations. This will be done inside a loop:

```r
for (Each.Mutation in as.character(Single.Mutations)) {
    # some code here
}
```
I'll go through the code inside the loop next. We first extract the Final PSI vs Starting PSI dataframe from `Final.Vs.Starting.PSI.List` and filter it so it only includes low-variance genotypes as well as genotypes whose inclusion is predicted to be less than or equal to 100%:

```r
# Extract the DF with effects from this mutation
Mutant.DF <- Final.Vs.Starting.PSI.List[[Each.Mutation]]

# Use low-variance genotypes for curve fitting
Low.Noise.Rows <- which(as.character(Mutant.DF$Genotype.Final) %in% Low.SD.Variants & as.character(Mutant.DF$Genotype.Starting) %in% Low.SD.Variants)
Mutant.DF <- Mutant.DF[Low.Noise.Rows,]

# and discard genotypes with predicted PSI > 100%
Impossible.Rows <- which(Mutant.DF$Final.PSI > 100 | Mutant.DF$Starting.PSI > 100)
Mutant.DF <- Mutant.DF[-Impossible.Rows,]
```
If the mutation only has one behaviour, then we generate its principal curve:

```r
# Principal curve for this mutation's dataset
PR.Curve <- prcurve(X = data.frame(X = Mutant.DF$Starting.PSI,
                                   Y = Mutant.DF$Final.PSI),
                    complexity = 4,
                    method = "pca")
```
And save it in the `Principal.Curves` list from before:

```r
# save principal curve
Principal.Curves[[Each.Mutation]] <- PR.Curve
```

But if the mutation has two behaviours, we first split `Mutant.DF` in two, depending on whether or not the epistatic mutation is present in the background:

```r
# which is the epistatic mutation?
Epistatic.Mutation <- as.character(Epistasis.Table$Epistatic)[which(as.character(Epistasis.Table$Main) == Each.Mutation)]

# rows with epistatic mutation
Epistatic.Rows <- which(sapply(as.character(Mutant.DF$Genotype.Starting),
                               function(x){
                                 Epistatic.Mutation %in% strsplit(x,";")[[1]]
                               }))

# subset of Mutant.DF only with rows containing epistatic mutation
Mutant.DF.1 <- Mutant.DF[Epistatic.Rows,]

# subset of Mutant.DF without rows containing epistatic mutation
Mutant.DF.2 <- Mutant.DF[-Epistatic.Rows,]
```
A principal curve is generated for each of these two data frames:

```r
# Principal curve for 1st behaviour
PR.Curve.1 <- prcurve(X = data.frame(X = Mutant.DF.1$Starting.PSI,
                                     Y = Mutant.DF.1$Final.PSI),
                      complexity = 4,
                      method = "pca")
                      
# Principal curve for 2nd behaviour
PR.Curve.2 <- prcurve(X = data.frame(X = Mutant.DF.2$Starting.PSI,
                                     Y = Mutant.DF.2$Final.PSI),
                      complexity = 4,
                      method = "pca")
```
The two principal curves are then saved in the `Principal.Curves` list:

```r
# save principal curve 1
This.Mutation.Name <- paste(Each.Mutation,
                            ".1",
                            sep = "",
                            collapse = "")
Principal.Curves[[This.Mutation.Name]] <- PR.Curve.1

# save principal curve 2
This.Mutation.Name <- paste(Each.Mutation,
                            ".2",
                            sep = "",
                            collapse = "")
Principal.Curves[[This.Mutation.Name]] <- PR.Curve.2
```
After all the principal curves have been generated, we will use them to estimate each mutation's maximum effect size and the starting PSI at which that effect size is observed:

```r
Max.Starting <- c()
Max.Deltas <- c()

for (each.prcurve in Principal.Curves){
  Index <- which(abs(each.prcurve$s[each.prcurve$tag,][,2] - each.prcurve$s[each.prcurve$tag,][,1]) == max(abs(each.prcurve$s[each.prcurve$tag,][,2] - each.prcurve$s[each.prcurve$tag,][,1])))[1]
  Max.Starting <- c(Max.Starting, each.prcurve$s[each.prcurve$tag,][Index,1])
  Max.Deltas <- c(Max.Deltas, each.prcurve$s[each.prcurve$tag,][Index,2] - each.prcurve$s[each.prcurve$tag,][Index,1])
}
```
We can plot the results:

```r
plot(NULL,
     xlim = c(-100,100),
     ylim = c(0,100),
     xlab = "",
     ylab = "",
     axes = F)
abline(50,-0.5,
       lwd = 2,
       lty = 2,
       col = "gray80")
par(new=T)
plot(y=Max.Starting,
     x=Max.Deltas,
     ylim = c(0,100),
     xlim = c(-100,100),
     pch = 19,
     col =  rgb(0,0,0,0.7),
     cex = 1.5,
     las = 1,
     xlab = "Max effect (Delta PSI)",
     ylab = "Starting PSI")
```
<p align="center">
  <img width = 450 height = 450 src="Figures/004_startingPSI_vs_max_effect2.png">
  <br> Figure 4F
</p>

## 5. Effect of tau

## 1. Assuming tau = 0
