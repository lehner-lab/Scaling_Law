#########################
## 1. Data preparation ##
#########################


# load the file containing the entire dataset
load("002_whole_dataset.RData")

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

# start an empty list of length 12
Final.Vs.Starting.PSI.List <- vector(mode = "list",
                                     length = length(Single.Mutations))

# each element in the list is named after one of the 12 mutations
names(Final.Vs.Starting.PSI.List) <- as.character(Single.Mutations)

# loop over the 12 mutations
for (i in 1:(length(Single.Mutations))) {
  # which is this mutation?
  My.ID <- Single.Mutations[i]
  
  # keep track of the iteration we're in...
  print(My.ID)
  
  # Get the rows containing my target mutation
  Rows.With.My.ID.Plus.SomethingElse <- sapply(X = as.character(Whole.Dataset$Mutation.IDs),
                                               FUN = function(x){
                                                 x.vector <- strsplit(x, ";")[[1]]
                                                 My.ID %in% x.vector
                                               })
  
  # get final psis
  Final.PSI <- Whole.Dataset$Mean[Rows.With.My.ID.Plus.SomethingElse]
  Low.Final.PSI <- Whole.Dataset$Low.PSI[Rows.With.My.ID.Plus.SomethingElse]
  High.Final.PSI <- Whole.Dataset$High.PSI[Rows.With.My.ID.Plus.SomethingElse]
  
  # Extract all those "and something else" that are in the background
  Mutations.In.Background <- sapply(X = as.character(Whole.Dataset$Mutation.IDs)[Rows.With.My.ID.Plus.SomethingElse],
                                    FUN = function(x){
                                      x.vector <- strsplit(x, ";")[[1]]
                                      which.is.my.id <- which(x.vector == My.ID)
                                      x.vector <- x.vector[-which.is.my.id]
                                      if (length(x.vector) == 0){
                                        x.vector <- ""
                                      }
                                      x.vector
                                    })
  
  
  # Find the rows that contain ONLY those "something else" IDs
  Index.Starting.Sequence <- lapply(X = Mutations.In.Background,
                                    FUN = function(x){
                                      x <- as.character(x)
                                      ID.To.Search <- paste(x, collapse = ";")[[1]]
                                      
                                      which(as.character(Whole.Dataset$Mutation.IDs) == ID.To.Search)
                                    })
  # unlist
  Index.Starting.Sequence <- unlist(Index.Starting.Sequence)
  
  # get starting psis
  Starting.PSI <- Whole.Dataset$Mean[Index.Starting.Sequence]
  Low.Starting.PSI <- Whole.Dataset$Low.PSI[Index.Starting.Sequence]
  High.Starting.PSI <- Whole.Dataset$High.PSI[Index.Starting.Sequence]
  
  
  # delta psi is final minus starting psi
  Delta.PSI <- Final.PSI - Starting.PSI
  
  
  #save as a data frame
  This.ID.Dataframe <- data.frame(Genotype.Final = as.character(Whole.Dataset$Mutation.IDs)[Rows.With.My.ID.Plus.SomethingElse],
                                  Genotype.Starting = as.character(Whole.Dataset$Mutation.IDs)[Index.Starting.Sequence],
                                  Final.PSI = Final.PSI,
                                  Low.Final.PSI = Low.Final.PSI,
                                  High.Final.PSI = High.Final.PSI,
                                  Starting.PSI = Starting.PSI,
                                  Low.Starting.PSI = Low.Starting.PSI,
                                  High.Starting.PSI = High.Starting.PSI,
                                  Delta.PSI = Delta.PSI)
  
  Final.Vs.Starting.PSI.List[[i]] <- This.ID.Dataframe
}

# save the file
save(Final.Vs.Starting.PSI.List, file = "004_final_vs_starting_psi_list.RData")




#################################################
## 2. Background where a mutation can be added ##
#################################################

# load files
load("002_whole_dataset.RData")
load("004_final_vs_starting_psi_list.RData")


# T19G is in the 4th position in the list
This.DF <- Final.Vs.Starting.PSI.List[[4]]

# change row names from sequence to genotype id
rownames(Whole.Dataset) <- as.character(Whole.Dataset$Mutation.IDs)

# number of mutations in starting genotype id
Number.Of.Mutations <- Whole.Dataset[as.character(This.DF$Genotype.Starting), "No.Mutations"]


# standard deviation of PSI in starting genotypes
SD <- Whole.Dataset[as.character(This.DF$Genotype.Starting), "SD"]


# table for bar plot
Table.Barplot <- rbind(
  c(1, as.vector(table(Number.Of.Mutations[which(SD < 10)])), 0),
  c(1, as.vector(table(Number.Of.Mutations))) - c(1, as.vector(table(Number.Of.Mutations[which(SD < 10)])), 0)
)

###############################################################
# PLOT
#

# draw the bar plot
par(pty="s")
barplot(Table.Barplot,
        col = c("indianred1", "gray90"),
        border = NA,
        ylim = c(0,400),
        names.arg = 0:10,
        las = 1,
        ylab = "Frequency")
        
#
#
###############################################################





#########################################
## 3. Distribution of mutation effects ##
#########################################


# load files
load("002_whole_dataset.RData")
load("004_final_vs_starting_psi_list.RData")

# mutations in the order in which they occur along the exon
Names.In.Sequence.Order <- c("C-18-T",
                             "C-18-G",
                             "T-19-G",
                             "T-24-C",
                             "G-26-T",
                             "C-32-T",
                             "G-35-T",
                             "C-39-T",
                             "C-41-G",
                             "G-44-A",
                             "T-49-C",
                             "G-51-C"
)

# before plotting, get the IDs of genotypes with low variance
Low.Noise.Genotypes <- as.character(Whole.Dataset$Mutation.IDs)[which(Whole.Dataset$SD < 10)]

###############################################################
# PLOT
#

par(mar = c(0,0,1,1) + 0.1)
par(oma = c(5,4,0,0) + 0.1)
par(mfrow=c(3,4))
par(pty="s")

# loop through the 12 mutations
for (i in 1:12){
  This.Number <- i
  
  # what is the id of this mutation?
  This.Mutation <- as.character(Names.In.Sequence.Order)[This.Number]
  
  # extract the data frame
  This.Mutation.DF <- Final.Vs.Starting.PSI.List[[This.Mutation]]
  
  # Get low noise rows
  Low.Noise.Rows <- which(as.character(This.Mutation.DF$Genotype.Final) %in% Low.Noise.Genotypes & as.character(This.Mutation.DF$Genotype.Starting) %in% Low.Noise.Genotypes)
  
  # subset This.Mutation.DF
  This.Mutation.DF <- This.Mutation.DF[Low.Noise.Rows,]
  
  # draw the histogram
  if (! This.Mutation %in% c("C-32-T", "T-24-C", "C-41-G")) {
    hist(This.Mutation.DF$Delta.PSI,
         xlim = c(-100,100),
         ylim = c(0,100),
         col = "gray80",
         border = "gray80",
         ann = F,
         axes = F,
         breaks = 10)
  } else {
    hist(This.Mutation.DF$Delta.PSI,
         xlim = c(-100,100),
         ylim = c(0,100),
         col = "gray80",
         border = "gray80",
         ann = F,
         axes = F,
         breaks = 5)
  }
  
  # mutation ID
  text(x = -100,
       y = 90,
       labels = paste(strsplit(This.Mutation, split = "-")[[1]], sep = "", collapse = ""),
       pos = 4, cex = 1)
  
  # surround the plot with a box (no axes were drawn yet)
  box()
  
  # draw left-hand side axis if this plot is to the left
  if (i %in% c(1,5,9)){
    axis(side = 2, las = 1)
  }
  
  # draw bottom axis if this plot is at the bottom
  if (i %in% c(9,10,11,12)) {
    axis(side = 1)
  }
  
}

# axis labels
title(xlab = "Delta PSI",
      ylab = "Frequency",
      outer = TRUE, line = 3)


par(mar = c(5,4,4,2)+0.1)
par(oma = c(0,0,0,0))
par(mfrow=c(1,1))

#
#
###############################################################



###############################################################
# PLOT
#

par(mar = c(0,0,1,1) + 0.1)
par(oma = c(5,4,0,0) + 0.1)
par(mfrow=c(3,4))
par(pty="s")

for (i in 1:12){
  This.Number <- i
  
  # what is the id of this mutation?
  This.Mutation <- as.character(Names.In.Sequence.Order)[This.Number]
  
  # extract the data frame
  This.Mutation.DF <- Final.Vs.Starting.PSI.List[[This.Mutation]]
  
  # draw the histogram
  if (! This.Mutation %in% c("C-32-T", "T-24-C", "C-41-G", "G-35-T", "G-44-A")) {
    hist(This.Mutation.DF$Delta.PSI,
         xlim = c(-100,100),
         ylim = c(0,600),
         col = "gray80",
         border = "gray80",
         ann = F,
         axes = F,
         breaks = 10)
  } else {
    hist(This.Mutation.DF$Delta.PSI,
         xlim = c(-100,100),
         ylim = c(0,600),
         col = "gray80",
         border = "gray80",
         ann = F,
         axes = F,
         breaks = 20)
  }
  
  # write mutation ID
  text(x = -100,
       y = 540,
       labels = paste(strsplit(This.Mutation, split = "-")[[1]], sep = "", collapse = ""),
       pos = 4, cex = 1)
  
  # surround the plot with a box (no axes were drawn yet)
  box()
  
  # draw left-hand side axis if this plot is to the left
  if (i %in% c(1,5,9)){
    axis(side = 2, las = 1)
  }
  
  # draw bottom axis if this plot is at the bottom
  if (i %in% c(9,10,11,12)) {
    axis(side = 1)
  }
  
}

# axis labels
title(xlab = "Delta PSI",
      ylab = "Frequency",
      outer = TRUE, line = 3)


par(mar = c(5,4,4,2)+0.1)
par(oma = c(0,0,0,0))
par(mfrow=c(1,1))

#
#
###############################################################

