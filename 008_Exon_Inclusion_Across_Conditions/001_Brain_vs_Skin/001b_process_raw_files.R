# load library
library(psichomics)

# Load Human (hg19/GRCh37 assembly) annotation
human <- listSplicingAnnotations()[[1]]
annotation <- loadAnnotation(human)

# Events to quantify (just look at alt exons for now)
eventType <- getSplicingEventTypes()[1]

# minimum read counts
minReads <- 10

# where is my file?
Junction.File.Path <- "GTEx_Analysis_2016-01-15_v7_STARv2.4.2a_junctions.gct"

# load it
GTEx.Junctions <- loadGtexData(junctionQuant = Junction.File.Path)
junctionQuant <- GTEx.Junctions$GTEx$"Junction quantification"

# Quantify
PSI.Estimates <- quantifySplicing(annotation = annotation,
                                  junctionQuant = junctionQuant,
                                  eventType = eventType, 
                                  minReads = minReads)

# Load table with sample information
GTEX.Dataset <- read.table(file = "IDs_Tissues.txt",
                           header = T,
                           sep = "\t")
# rename columns & rows
colnames(GTEX.Dataset) <- c("ID", "Tissue")
rownames(GTEX.Dataset) <- as.character(GTEX.Dataset$ID)

# find samples that are found in the annotations data frame and
# the PSI estimates dataframe
Common.Samples <- intersect(as.character(rownames(GTEX.Dataset)),
                            as.character(colnames(PSI.Estimates)))

# subset
GTEX.Dataset <- GTEX.Dataset[Common.Samples,]
PSI.Estimates <- PSI.Estimates[,Common.Samples]

save(GTEX.Dataset, PSI.Estimates, file = "Tissues_Compared_Datasets.RData")