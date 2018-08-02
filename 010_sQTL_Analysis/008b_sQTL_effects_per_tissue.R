####################
## 1. Preparation ##
####################

# load data handling lib
library(data.table)

# sQTL ID is fed from the command line
args = commandArgs(trailingOnly=TRUE)
sQTL.ID <- args[1]

# splicing event ID taken from the BED file
Splicing.Event.ID <- as.character(read.table("Gene_Regions.bed")[1,4])

# load GTEX.Dataset.Annotations and PSI.Estimates
load("../../Data/Annotations_PSIestimates.RData")

# import a vector with the IDs of all 35153 splicing events analysed
All.Splicing.Events <- as.character(read.table("../../Gene_Regions.bed", header = T)$ID)

# which row contains our splicing event of interest?
Splicing.Event.Row <- which(All.Splicing.Events == Splicing.Event.ID)

# PSI values for this splicing event across all our samples
This.Event.PSI <- as.numeric(PSI.Estimates[Splicing.Event.Row,])

# Build a table with information about this splicing event
This.Splicing.Event.Table <- GTEX.Dataset.Annotations
This.Splicing.Event.Table$PSI <- This.Event.PSI
This.Splicing.Event.Table$Subject.ID <- gsub(pattern = "\\.",
                                             replacement = "-",
                                             x = This.Splicing.Event.Table$Subject.ID)

# load parsed genotypes file
Genotypes <- fread(input = "genotypes.MAF01.GeneRegions.Parsed.txt",
                   skip = "#CHROM")

# row for this sQTL
sQTL.Row <- which(as.character(Genotypes$ID) == sQTL.ID)

# extract genotypes for this sQTL
This.Splicing.Event.Genotypes <- as.numeric(Genotypes[sQTL.Row,4:ncol(Genotypes)])
names(This.Splicing.Event.Genotypes) <- colnames(Genotypes)[4:ncol(Genotypes)]

# empty vector which we'll fill in with genotype data, in
# the same order as the samples from This.Splicing.Event.Table
Genotypes.For.Table <- c()

# fill in the empty vector with sQTL genotype in the correct order
for(each.id in This.Splicing.Event.Table$Subject.ID) {
  if (each.id %in% names(This.Splicing.Event.Genotypes)) {
    Genotypes.For.Table <- c(Genotypes.For.Table, as.numeric(This.Splicing.Event.Genotypes[each.id]))
  } else {
    Genotypes.For.Table <- c(Genotypes.For.Table, NA)
  }
}

This.Splicing.Event.Table$Genotype <- Genotypes.For.Table




###############################
## 2. Calculate sQTl effects ##
###############################

# get all possible tissues
All.Tissues <- unique(GTEX.Dataset.Annotations$Tissue)

# initialise empty vector for PSI at genotype 0
Vector.Tissue.PSI.At.Gt.0 <- rep(NA, length(All.Tissues))
names(Vector.Tissue.PSI.At.Gt.0) <- as.character(All.Tissues)

# empty vector for sQTL effect
Vector.Tissue.Effect.Of.QTL <- rep(NA, length(All.Tissues))
names(Vector.Tissue.Effect.Of.QTL) <- as.character(All.Tissues)

# loop through the different tissues and get:
# - PSI of splicing event at gt == 0
# - Effect of sQTL per unit gt in all tissues
for (each.tissue in All.Tissues) {
  SubDF <- This.Splicing.Event.Table[which(This.Splicing.Event.Table$Tissue == each.tissue),]
  
  Which.Genotypes.To.Consider <- which(table(SubDF$Genotype) >= 7)
  
  if (length(Which.Genotypes.To.Consider) >= 2) {
    Genotype.To.Ignore <- which(table(SubDF$Genotype) < 7) - 1
    
    if (! length(Genotype.To.Ignore) == 0){
      SubDF <- SubDF[- which(SubDF$Genotype == Genotype.To.Ignore),]
    }
    
    PSI.At.Gt.0 <- lm(PSI ~ Genotype, data = SubDF)$coefficients[1]
    Vector.Tissue.PSI.At.Gt.0[each.tissue] <- PSI.At.Gt.0
    
    Effect.Of.QTL <- lm(PSI ~ Genotype, data = SubDF)$coefficients[2]
    Vector.Tissue.Effect.Of.QTL[each.tissue] <- Effect.Of.QTL
  }
  
}

# join everything into one long vector
This.sQTL.Information <- c(Vector.Tissue.PSI.At.Gt.0,
                           Vector.Tissue.Effect.Of.QTL)

# save the vector
save(This.sQTL.Information, file = paste(paste("Effect.Of.sQTL.In.Different.Tissues",
                                               sQTL.ID,
                                               sep = "@"),
                                         ".RData",
                                         sep = ""))
