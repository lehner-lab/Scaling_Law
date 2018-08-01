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

# load vector with allowed subjects (common to genotypes + junctions file)
Allowed.Subjects <- as.character(read.table(file = "Common_Subjects.txt")$V1)

# all Samples from the PSI data frame
All.Samples <- colnames(PSI.Estimates)

# which of these samples can I work with?
Can.I.Work.With.This.Sample <- sapply(All.Samples,
                                      function(x){
                                        this.subject <- paste(strsplit(x = x,
                                                                       split = "\\.")[[1]][1:2],
                                                              collapse = "-")
                                        this.subject %in% Allowed.Subjects
                                      })

Allowed.Samples <- which(Can.I.Work.With.This.Sample)

# update PSI.Estimates so it only contains samples we can work with
PSI.Estimates <- PSI.Estimates[,Allowed.Samples]

# Look at mean NA value per splicing event
NAs.Per.Splicing.Event <- apply(X = PSI.Estimates,
                                MARGIN = 1,
                                FUN = function(x) mean(is.na(x)))

# select only events where I don't have any NAs
Not.All.NAs <- which(NAs.Per.Splicing.Event < 1)
PSI.Estimates <- PSI.Estimates[Not.All.NAs,]

# save
save(PSI.Estimates, file = "PSI.Estimates.RData")
