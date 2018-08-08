
##############################
## 1. Splicing events table ##
##############################

# Load massive data frame with all the psi values
load("PSI.Estimates.RData")

# chromosome id
Chromosomes <- sapply(as.character(rownames(PSI.Estimates)),
                      function(x){
                        strsplit(x,"_")[[1]][2]
                      })

# gene symbol
Genes <- sapply(as.character(rownames(PSI.Estimates)),
                function(x){
                  split_id <- strsplit(x,"_")[[1]]
                  split_id[length(split_id)]
                })

# build a table with the following information:
# - chromosome (column 1)
# - id of splicing event (column 2)
# - gene symbol where splicing event is found
Splicing.Events.Table <- data.frame(Chr = Chromosomes,
                                    ID = as.character(rownames(PSI.Estimates)),
                                    Gene = Genes)








########################
## 2. Build BED table ##
########################


# some data handling libraries
library(data.table)
library(dplyr)


# load Gencode annotations file using 'fread' function
Gencode.Annotations <- fread(input = "table.gencode.v19.annotation.gff3")

# add column names
colnames(Gencode.Annotations) <- c("Chromosome",
                                   "Source",
                                   "Type",
                                   "Start",
                                   "End",
                                   "Score",
                                   "Strand",
                                   "Phase",
                                   "Attributes")

# we're only interested in genes
Gencode.Annotations.Genes <- Gencode.Annotations[Type == "gene"]

# extract gene symbols
Gencode.Annotations.Genes$Gene <- sapply(Gencode.Annotations.Genes$Attributes,
                                         function(x){
                                           strsplit(strsplit(x,";")[[1]][6], "=")[[1]][2]
                                         })


# add gene start positions to Splicing.Events.Table
Splicing.Events.Table$Start <- sapply(Splicing.Events.Table$Gene,
                                function(x){
                                  
                                  row.number <- which(Gencode.Annotations.Genes$Gene == x)
                                  if (length(row.number) > 1) {
                                    row.number <- row.number[1]
                                  }
                                  
                                  if (! length(row.number) == 0){
                                    Start <- Gencode.Annotations.Genes$Start[row.number] - 1
                                  } else {
                                    Start <- NA
                                  }
                                  
                                  Start
                                  
                                })

# add gene end positions to Splicing.Events.Table
Splicing.Events.Table$End <- sapply(Splicing.Events.Table$Gene,
                              function(x){
                                
                                row.number <- which(Gencode.Annotations.Genes$Gene == x)
                                if (length(row.number) > 1) {
                                  row.number <- row.number[1]
                                }
                                
                                if (! length(row.number) == 0){
                                  End <- Gencode.Annotations.Genes$End[row.number]
                                } else {
                                  End <- NA
                                }
                                
                                End
                                
                              })


# shuffle Splicing.Events.Table columns around to fit BED format
Splicing.Events.Table <- Splicing.Events.Table %>% select(Chr, Start, End, ID, Gene)
Splicing.Events.Table <- Splicing.Events.Table[complete.cases(Splicing.Events.Table),]


# save BED file
write.table(x = Splicing.Events.Table,
            file = "Gene_Regions.bed",
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")
