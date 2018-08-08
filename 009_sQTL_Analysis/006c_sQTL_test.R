
#################################
## 1. Load libraries and files ##
#################################

# library to perform likelihood-ratio test
library(lmtest)

# data table library for data handling
library(data.table)

# Load PSI.Estimates and GTEx annotations
load("../../Filtered_Annotations_PSIestimates.RData")

# load bed file with information about this particular splicing event
Splicing.Events <- fread(input = "Gene_Regions.bed")
colnames(Splicing.Events) <- c("Chr", "Start", "End", "ID", "Gene")

# load parsed genotype file for this particular splicing event
Genotypes <- fread("genotypes.MAF01.GeneRegions.Parsed.txt",
                   skip = "#CHROM")
colnames(Genotypes) <- sapply(X = colnames(Genotypes),
                              FUN = function(x){
                                gsub(pattern = "-", replacement = "\\.", x = x)
                              })

# filter genotypes file to keep only subjects included in our analysis
Subjects.We.Can.Look.At <- which(colnames(Genotypes) %in% GTEX.Dataset.Annotations$Subject.ID)
Genotypes <- Genotypes[,.SD,,.SDcols=c(1,2,3,Subjects.We.Can.Look.At)]




######################
## 2. Test for sQTL ##
######################

# inclusion levels for this particular exon skipping event
This.Splicing.Event <- as.character(Splicing.Events$ID)
This.Event.PSI <- as.numeric(PSI.Estimates[This.Splicing.Event,])

# loop through all the potential sQTLs and perform a statistical test if appropriate
if (nrow(Genotypes) != 0) {
  # build a df for this particular exon splicing event
  This.Splicing.Event.Table <- GTEX.Dataset.Annotations
  This.Splicing.Event.Table$PSI <- This.Event.PSI
  
  # start an empty list
  sQTL.List <- vector(mode = "list", length = nrow(Genotypes))
  
  # loop through all the potential QTLs and perform a statistical test if appropriate
  for (j in 1:nrow(Genotypes)){
    
    # start an empty vector
    This.Splicing.Event.Genotypes <- c()
    
    # fill in empty vector with genotypes (0/1/2)
    for(each.id in This.Splicing.Event.Table$Subject.ID) {
      if (each.id %in% colnames(Genotypes)) {
        This.Splicing.Event.Genotypes <- c(This.Splicing.Event.Genotypes,
                                           as.numeric(Genotypes[j,
                                                                each.id,
                                                                with=FALSE]))
      } else {
        This.Splicing.Event.Genotypes <- c(This.Splicing.Event.Genotypes,
                                           NA)
      }
    }
    
    # put genotypes in new column of This.Splicing.Event.Table
    This.Splicing.Event.Table$Genotype <- This.Splicing.Event.Genotypes
    
    # what tissues do I have?
    Vector.of.Tissues.Analysed <- unique(as.character(This.Splicing.Event.Table$Tissue))
    
    # start an empty list for p values sorted by tissue
    Tissues.List <- vector(mode = "list",
                           length = length(Vector.of.Tissues.Analysed))
    
    # elements of the list named after the different tissues in our dataset
    names(Tissues.List) <- Vector.of.Tissues.Analysed
    
    for (each.tissue in Vector.of.Tissues.Analysed) {
      
      # a sub-dataframe only with rows corresponding to this tissue
      SubDF <- This.Splicing.Event.Table[which(as.character(This.Splicing.Event.Table$Tissue) == each.tissue),]
      SubDF <- SubDF[complete.cases(SubDF),]
      
      # only do a statistical test if we have at least two of (0, 1, 2);
      # and if the no. of observations in at least two of the categories is > 7
      if (length(table(SubDF$Genotype)) >= 2 & sum(table(SubDF$Genotype) > 7) >= 2) {
        
        # value to be predicted by the model
        lm.psi <- SubDF$PSI
        lm.psi <- round(lm.psi*100)
        lm.psi <- cbind(lm.psi, 100-lm.psi)
        
        # genotype predictor
        lm.genotype <- as.numeric(as.character(SubDF$Genotype))
        
        # ischaemic time predictor
        # (divide by 1000 because function crashes with large numbers)
        lm.ischaemic <- SubDF$Ischaemic.Time/1000
        
        # sex predictor
        lm.sex <- as.numeric(as.character(SubDF$Sex))
        
        # age predictor
        # (divide by 100 since function crashes with large numbers)
        lm.age <- as.numeric(as.character(SubDF$Age))/100
        
        # model including genotype variable
        res<-try(suppressMessages(glm(lm.psi ~ lm.genotype + lm.ischaemic + lm.sex + lm.age,
                                      family=binomial)),
                 silent=TRUE)
        
        # model without genotype variable
        res0<-try(suppressMessages(glm(lm.psi ~ lm.ischaemic + lm.sex + lm.age,
                                       family=binomial)),
                  silent=TRUE)
        
        # do a likelihood-ratio test and take p value
        P.Value <- lrtest(res,res0)$"Pr(>Chisq)"[2]
        Estimate <- summary(res)$coefficients[2,1]
        
        Tissues.List[[each.tissue]] <- c(P.Value, Estimate)
      } else {
        Tissues.List[[each.tissue]] <- NA
      }
    }
    
    # put Tissues.List inside sQTL.List
    sQTL.List[[j]] <- Tissues.List
    
  }
  
  # name elements of sQTL.List with the sQTL IDs
  names(sQTL.List) <- Genotypes$ID
  
  # save sQTL.List
  save(sQTL.List, file = "sQTL.RData")
  
}
