
##############################################
## 1. Genotypes different from the ancestor ##
##############################################

# load Whole.Dataset
load("002_whole_dataset.RData")

# low variance subset of the data
Low.Variance <- Whole.Dataset[which(Whole.Dataset$SD < 10),]

# wilcoxon test for difference with 96% PSI
P.Different.From.Ancestor <- apply(X = Low.Variance[,2:10],
                                   MARGIN = 1,
                                   FUN = function(x){
                                     wilcox.test(x = x, mu = 96)$p.value
                                   })

# fdr-correct the p values
FDR.Different.From.Ancestor <- p.adjust(p = P.Different.From.Ancestor,
                                        method = "fdr")

# count how many are different from ancestor
Number.Significant.Differences <- sum(FDR.Different.From.Ancestor < 0.05)
Total.Genotypes <- length(FDR.Different.From.Ancestor)

Number.Significant.Differences/Total.Genotypes * 100




####################################################
## 2. Genotypes different from the human sequence ##
####################################################

# wilcoxon test for difference with 59.5% PSI
P.Different.From.Human <- apply(X = Low.Variance[,2:10],
                                MARGIN = 1,
                                FUN = function(x){
                                  wilcox.test(x = x, mu = 59.5)$p.value
                                })

# fdr-correct the p values
FDR.Different.From.Human <- p.adjust(p = P.Different.From.Human,
                                     method = "fdr")

# count how many are different from human
Number.Significant.Differences <- sum(FDR.Different.From.Human < 0.05)
Total.Genotypes <- length(FDR.Different.From.Human)

Number.Significant.Differences/Total.Genotypes * 100




##########################
## 3. Behaviour of T19G ##
##########################

# load final PSI - starting PSI dataset
load("004_final_vs_starting_psi_list.RData")

# extract table for mutation T-19-G
T19G.Table <- Final.Vs.Starting.PSI.List[["T-19-G"]]

# factor to character
T19G.Table$Genotype.Final <- as.character(T19G.Table$Genotype.Final)
T19G.Table$Genotype.Starting <- as.character(T19G.Table$Genotype.Starting)

# change the ID of the ancestral genotype
T19G.Table$Genotype.Starting[which(T19G.Table$Genotype.Starting == "")] <- "ANCESTOR"

# set Whole.Dataset row names to be the genotype ID
rownames(Whole.Dataset) <- as.character(Whole.Dataset$Mutation.IDs)
rownames(Whole.Dataset)[which(Whole.Dataset$Mutation.IDs == "")] <- "ANCESTOR"

# take all 9 replicates of each starting and final PSI in T19G.Table
Starting.PSI <- Whole.Dataset[as.character(T19G.Table$Genotype.Starting),2:10]
Final.PSI <- Whole.Dataset[as.character(T19G.Table$Genotype.Final),2:10]

# make a data frame with the difference
Delta.PSI <- Final.PSI - Starting.PSI

# calculate P values
P.Values.Delta.PSI <- apply(X = Delta.PSI,
                            MARGIN = 1,
                            FUN = function(x){
                              wilcox.test(x = as.numeric(x), mu = 0)$p.value
                            })

# FDR-correct
FDR.Values.Delta.PSI <- p.adjust(p = P.Values.Delta.PSI, method = "fdr")

# how many backgrounds promote inclusion/skipping?
sum(FDR.Values.Delta.PSI < 0.05 & T19G.Table$Delta.PSI > 0)
sum(FDR.Values.Delta.PSI < 0.05 & T19G.Table$Delta.PSI < 0)
