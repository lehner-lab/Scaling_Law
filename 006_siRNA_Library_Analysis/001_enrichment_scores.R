
##################
## 1. Load data ##
##################

# load siControl replicates
for (i in 1:3) {
  # the name of the file I need to read
  This.File.Name <- paste("Control_RNAi_Rep_", i, ".counts", sep = "")
  # load it in R with a generic variable name
  This.Data.Frame <- read.table(This.File.Name, header = T)
  # set column names
  colnames(This.Data.Frame) <- c("Sequence", "Counts", "Mutations")
  # set row names
  rownames(This.Data.Frame) <- as.character(This.Data.Frame$Sequence)
  # variable name I want to assign it to
  This.Variable.Name <- paste("siControl", i, sep = "")
  # give it the new variable name
  assign(x = This.Variable.Name, value = This.Data.Frame)
}

# load siSF3B replicates
for (i in 1:3) {
  # the name of the file I need to read
  This.File.Name <- paste("SF3B_RNAi_Rep_", i, ".counts", sep = "")
  # load it in R with a generic variable name
  This.Data.Frame <- read.table(This.File.Name, header = T)
  # set column names
  colnames(This.Data.Frame) <- c("Sequence", "Counts", "Mutations")
  # set row names
  rownames(This.Data.Frame) <- as.character(This.Data.Frame$Sequence)
  # variable name I want to assign it to
  This.Variable.Name <- paste("siSF3B", i, sep = "")
  # give it the new variable name
  assign(x = This.Variable.Name, value = This.Data.Frame)
}


# load input
for (i in 1:3) {
  # the name of the file I need to read
  This.File.Name <- paste("Doped_Library_Input_Rep_", i, ".counts", sep = "")
  # load it in R with a generic variable name
  This.Data.Frame <- read.table(This.File.Name, header = T)
  # set column names
  colnames(This.Data.Frame) <- c("Sequence", "Counts", "Mutations")
  # set row names
  rownames(This.Data.Frame) <- as.character(This.Data.Frame$Sequence)
  # variable name I want to assign it to
  This.Variable.Name <- paste("Input", i, sep = "")
  # give it the new variable name
  assign(x = This.Variable.Name, value = This.Data.Frame)
}

# find sequences common to all experiments
Common.Genotypes <- Reduce(f = intersect,
                           x = list(  as.character(siControl1$Sequence),
                                      as.character(siControl2$Sequence),
                                      as.character(siControl3$Sequence),
                                      as.character(siSF3B1$Sequence),
                                      as.character(siSF3B2$Sequence),
                                      as.character(siSF3B3$Sequence),
                                      as.character(Input1$Sequence),
                                      as.character(Input2$Sequence),
                                      as.character(Input3$Sequence)
                                      ))

# take only variants common to all datasets
siControl1 <- siControl1[Common.Genotypes,]
siControl2 <- siControl2[Common.Genotypes,]
siControl3 <- siControl3[Common.Genotypes,]

siSF3B1 <- siSF3B1[Common.Genotypes,]
siSF3B2 <- siSF3B2[Common.Genotypes,]
siSF3B3 <- siSF3B3[Common.Genotypes,]

Input1 <- Input1[Common.Genotypes,]
Input2 <- Input2[Common.Genotypes,]
Input3 <- Input3[Common.Genotypes,]

# cleanup
rm(This.Data.Frame,
   This.Variable.Name,
   This.File.Name,
   i)




####################################
## 2. Calculate enrichment scores ##
####################################

All.Expts <- data.frame(Input1 = Input1[rownames(Input1), "Counts"],
                        Input2 = Input2[rownames(Input1), "Counts"],
                        Input3 = Input3[rownames(Input1), "Counts"],
                        siControl1 = siControl1[rownames(Input1), "Counts"],
                        siControl2 = siControl2[rownames(Input1), "Counts"],
                        siControl3 = siControl3[rownames(Input1), "Counts"],
                        siSF3B1 = siSF3B1[rownames(Input1), "Counts"],
                        siSF3B2 = siSF3B2[rownames(Input1), "Counts"],
                        siSF3B3 = siSF3B3[rownames(Input1), "Counts"],
                        Sequence = as.character(Input1[rownames(Input1), "Sequence"]),
                        Mutations = Input1[rownames(Input1), "Mutations"])
All.Expts$Sequence <- as.character(All.Expts$Sequence)
rownames(All.Expts) <- rownames(Input1)

# Frequency of each sequence in the input
Input.Frequency <- apply(X = All.Expts[,1:3],
                         MARGIN = 2,
                         FUN = function(x){
                           return(x/sum(x, na.rm = T))
                         })

# Frequency of each sequence in the siControl
siControl.Frequency <- apply(X = All.Expts[,4:6],
                             MARGIN = 2,
                             FUN = function(x){
                               return(x/sum(x, na.rm = T))
                             })

# Frequency of each sequence in the siSF3B
siSF3B.Frequency <- apply(X = All.Expts[,7:9],
                          MARGIN = 2,
                          FUN = function(x){
                            return(x/sum(x, na.rm = T))
                          })

# Median frequency in the input
Input.Median.Frequency <- apply(X = Input.Frequency,
                                MARGIN = 1,
                                FUN = median, na.rm = T)

# calculate enrichment scores for siControl
Enrichment.Scores.siControl <- siControl.Frequency / Input.Median.Frequency
Enrichment.Scores.siControl <- as.data.frame(Enrichment.Scores.siControl)

# calculate enrichment scores for siSF3B
Enrichment.Scores.siSF3B <- siSF3B.Frequency / Input.Median.Frequency
Enrichment.Scores.siSF3B <- as.data.frame(Enrichment.Scores.siSF3B)

# save data frame
save(Enrichment.Scores.siControl,
     Enrichment.Scores.siSF3B,
     file = "001_enrichment_scores.RData")