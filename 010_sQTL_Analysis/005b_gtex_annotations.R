# data table package
library(data.table)

# load tables in R
All.Samples <- fread("IDs_Tissues_IschaemicTime.txt")
All.Subjects <- fread("Subject_Sex_Age.txt")

# subject ID column
All.Samples$SUBJID <- sapply(as.character(All.Samples$SAMPID),
                             function(x){
                               paste(strsplit(x, "\\.")[[1]][1:2],
                                     sep = ".",
                                     collapse = ".")
                             })

# note: for some reason no annotations for bone marrow subjects.
# Also ID has a different format. These samples are excluded in the merging step.
GTEX.Dataset.Annotations <- merge(All.Samples, All.Subjects)
colnames(GTEX.Dataset.Annotations) <- c("Subject.ID",
                                        "Sample.ID",
                                        "Tissue",
                                        "Ischaemic.Time",
                                        "Sex",
                                        "Age")

rownames(GTEX.Dataset.Annotations) <- as.character(GTEX.Dataset.Annotations$Sample.ID)


# load PSI.Estimates
load("PSI.Estimates.RData")

# What samples do we want to keep?
Sample.IDs.To.Keep <- intersect(as.character(colnames(PSI.Estimates)),
                                as.character(GTEX.Dataset.Annotations$Sample.ID))

# filter GTEX.Dataset.Annotations
GTEX.Dataset.Annotations <- GTEX.Dataset.Annotations[Sample.ID %in% Sample.IDs.To.Keep]

# make sure PSI.Estimates columns in the same order
PSI.Estimates <- PSI.Estimates[, as.character(GTEX.Dataset.Annotations$Sample.ID)]

# load splicing events table
Splicing.Events <- fread(input = "Gene_Regions.bed")

# make sure PSI.Estimates has the rows in the same order as Splicing.Events
PSI.Estimates <- PSI.Estimates[as.character(Splicing.Events$ID),]

# save
save(GTEX.Dataset.Annotations,
     PSI.Estimates,
     file="Filtered_Annotations_PSIestimates.RData")
