
##########################
## 1. Preparation stuff ##
##########################

# libraries i need
library(data.table)


# Load data
load("Annotations_PSIestimates.RData")
Vector.of.Tissues.Analysed <- unique(as.character(GTEX.Dataset.Annotations$Tissue))


# empty data table to fill with test results from previous step
Final.Data.Table <- data.table(Splicing.Event.ID = character(),
                               sQTL.ID = character(),
                               Tissue = character(),
                               P.Value = numeric())


for (i in 1:35153){
  print(i)

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
  
  # load BED file into R
  Splicing.Events <- fread(Gene.Regions.File)
  colnames(Splicing.Events) <- c("Chr", "Start", "End", "ID", "Gene")
  
  
  # Coordinates of gene associated to this splicing event
  This.Event.Start <- Splicing.Events$Start
  This.Event.End <- Splicing.Events$End
  
  # load genotypes file into R
  Genotypes <- fread(Genotypes.File,
                     skip = "#CHROM")
  
  # make sure snp positions fall within the gene boundaries.
  # because of how bed files are encoded, the last position should not be included
  Genotypes <- Genotypes[which(Genotypes$POS < This.Event.End & Genotypes$POS >= This.Event.Start), ]
  
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
  
  
  
  ################################
  ## 2. Gather sQTL information ##
  ################################
  
  # make sure the sQTL.RData exists; otherwise stop
  if (file.exists(File.Location)){
    
    # load sQTL.RData containing sQTL.List
    load(file = File.Location)
    
    # Subset sQTL to only include IDs found in Genotypes.
    # This step is necessary because I might have removed an ID
    # earlier when discussing the issues with BED format etc
    sQTL.List <- sQTL.List[Genotypes$ID]
    
    # make sure sQTL.List is not empty
    if (length(sQTL.List) > 0) {
      
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
      
      # build a data table with the four columns
      Temporary.Data.Table <- data.table(Splicing.Event.ID = Data.Frame.Splicing.ID,
                                         sQTL.ID = Data.Frame.QTL.ID,
                                         Tissue = Data.Frame.Tissue,
                                         P.Value = Data.Frame.P.Value)
      
      Final.Data.Table <- rbindlist(l = list(Final.Data.Table, Temporary.Data.Table))
    }
    
  }
  
}



################################################
## 3. Unique combination of significant sQTLs ##
################################################

# Bonferroni correction!!
Final.Data.Table <- Final.Data.Table[complete.cases(Final.Data.Table),]
Final.Data.Table$FDR.Value <- p.adjust(Final.Data.Table$P.Value, "bonferroni")

# subset of Final.Data.Table containing only significant rows
Significant.Final.Data.Table <- Final.Data.Table[which(Final.Data.Table$FDR.Value < 0.05),]

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

# save
write.table(x = Unique.Combos.DF,
            file = "Data/007_Significant_Splicing_Events_sQTLs_Combos.txt",
            col.names = FALSE,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)
