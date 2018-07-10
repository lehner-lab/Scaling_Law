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
par(pty="s")
plot(log2(Range.k6.To.Test),
     Calculated.PSIs*100,
     type = "l",
     col = "gray90",
     lwd = 3,
     xlab = "log2 k6",
     ylab = "PSI",
     las = 1)








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







# a bunch of maximum effect sizes
Max.Delta.PSIs <- seq(-100, 100, 20)

# the corresponding starting PSIs
Starting.PSIs <- 100 *(0.5 - 0.5 * Max.Delta.PSIs/100)

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











































######################
## principal curves ##
######################


# principal curves library
library(analogue)

# load datasets
load("../004_HEK293_Phylogeny_Library_Analysis/002_whole_dataset.RData")
load("../004_HEK293_Phylogeny_Library_Analysis/004_final_vs_starting_psi_list.RData")



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





# build a vector with low SD genotypes
Low.Noise.Dataset <- Whole.Dataset[which(Whole.Dataset$SD < 10),]
Low.SD.Variants <- as.character(Low.Noise.Dataset$Mutation.IDs)















Principal.Curves <- vector(mode = "list")




for (Each.Mutation in as.character(Single.Mutations)) {
  print(Each.Mutation)
  
  # Extract the DF with effects from this mutation
  Mutant.DF <- Final.Vs.Starting.PSI.List[[Each.Mutation]]
  
  # Use low-variance genotypes for curve fitting
  Low.Noise.Rows <- which(as.character(Mutant.DF$Genotype.Final) %in% Low.SD.Variants & as.character(Mutant.DF$Genotype.Starting) %in% Low.SD.Variants)
  Mutant.DF <- Mutant.DF[Low.Noise.Rows,]
  
  # and discard genotypes with predicted PSI > 100%
  Impossible.Rows <- which(Mutant.DF$Final.PSI > 100 | Mutant.DF$Starting.PSI > 100)
  Mutant.DF <- Mutant.DF[-Impossible.Rows,]
  
  
  # if 2 behaviours, build 2 curves
  if (Each.Mutation %in% Singles.With.Multiple.Paths){
    
    ## 1st behaviour!!!
    
    # which is the epistatic mutation?
    Epistatic.Mutation <- as.character(Epistasis.Table$Epistatic)[which(as.character(Epistasis.Table$Main) == Each.Mutation)]
    
    # rows with epistatic mutation
    Epistatic.Rows <- which(sapply(as.character(Mutant.DF$Genotype.Starting),
                                   function(x){
                                     Epistatic.Mutation %in% strsplit(x,";")[[1]]
                                   }))
    
    # subset of Mutant.DF only with rows containing epistatic mutation
    Mutant.DF.1 <- Mutant.DF[Epistatic.Rows,]
    
    # Principal curve for 1st behaviour
    PR.Curve.1 <- prcurve(X = data.frame(X = Mutant.DF.1$Starting.PSI,
                                         Y = Mutant.DF.1$Final.PSI),
                          complexity = 4,
                          method = "pca")
    
    # save principal curve 1
    This.Mutation.Name <- paste(Each.Mutation,
                                ".1",
                                sep = "",
                                collapse = "")
    Principal.Curves[[This.Mutation.Name]] <- PR.Curve.1
    
    
    
    ## 2nd behaviour!!
    
    # subset of Mutant.DF without rows containing epistatic mutation
    Mutant.DF.2 <- Mutant.DF[-Epistatic.Rows,]
    
    # Principal curve for 2nd behaviour
    PR.Curve.2 <- prcurve(X = data.frame(X = Mutant.DF.2$Starting.PSI,
                                         Y = Mutant.DF.2$Final.PSI),
                          complexity = 4,
                          method = "pca")
    
    # save principal curve 2
    This.Mutation.Name <- paste(Each.Mutation,
                                ".2",
                                sep = "",
                                collapse = "")
    Principal.Curves[[This.Mutation.Name]] <- PR.Curve.2
    
    
    
  } else { # if 1 behaviour, just one curve

    # Principal curve for this mutation's dataset
    PR.Curve <- prcurve(X = data.frame(X = Mutant.DF$Starting.PSI,
                                       Y = Mutant.DF$Final.PSI),
                        complexity = 4,
                        method = "pca")
    
    # save principal curve
    Principal.Curves[[Each.Mutation]] <- PR.Curve
  }
  
  
} # closing the Each.Mutation for loop






















Max.Starting <- c()
Max.Deltas <- c()

for (each.prcurve in Principal.Curves){
  
  Index <- which(abs(each.prcurve$s[each.prcurve$tag,][,2] - each.prcurve$s[each.prcurve$tag,][,1]) == max(abs(each.prcurve$s[each.prcurve$tag,][,2] - each.prcurve$s[each.prcurve$tag,][,1])))[1]
  Max.Starting <- c(Max.Starting, each.prcurve$s[each.prcurve$tag,][Index,1])
  Max.Deltas <- c(Max.Deltas, each.prcurve$s[each.prcurve$tag,][Index,2] - each.prcurve$s[each.prcurve$tag,][Index,1])
  
}

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





