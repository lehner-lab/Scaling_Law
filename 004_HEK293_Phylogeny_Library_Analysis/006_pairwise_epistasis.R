

# load Whole.Dataset
load("002_whole_dataset.RData")

# load Final.Vs.Starting.PSI
load("004_final_vs_starting_psi_list.RData")

# use Whole.Dataset to get the genotypes with low variance
Low.Noise.Genotypes <- as.character(Whole.Dataset$Mutation.IDs)[which(Whole.Dataset$SD < 10)]

# logit transform function
Logit.Transform <- function(x){
  log((x)/(1-x))
}


# initialise an empty data frame that I'll fill in with information
P.Values.Distance.DF <- data.frame(Main = c(),
                                   Epistatic = c(),
                                   Position.Main = c(),
                                   Position.Epistatic = c(),
                                   Distance = c(),
                                   Magnitude = c(),
                                   P.Value = c(),
                                   Minus.Log10.P.Value = c())





# loop through each mutations and test for interactions with all other mutations
for (k in 1:length(Final.Vs.Starting.PSI.List)){
  
  print(k)
  # Define Mutation
  This.Mutation <- names(Final.Vs.Starting.PSI.List)[k]
  
  # Make dataframe with before/after Mutation
  This.Mutation.DF <- Final.Vs.Starting.PSI.List[[This.Mutation]]
  
  # which rows in This.Mutation.DF contain only low-variance genotypes?
  Low.Noise.Indices <- which(as.character(This.Mutation.DF$Genotype.Final) %in% Low.Noise.Genotypes & as.character(This.Mutation.DF$Genotype.Starting) %in% Low.Noise.Genotypes)
  
  # filter high-variance genotypes out of This.Mutation.DF
  This.Mutation.DF <- This.Mutation.DF[Low.Noise.Indices,]
  
  # Find which other mutations the current mutation interacts with
  Other.Mutations <- names(Final.Vs.Starting.PSI.List)[which(names(Final.Vs.Starting.PSI.List) != This.Mutation)]
  
  # Cannot have 2 mutations in the same position
  if (This.Mutation == "C-18-G") {
    Other.Mutations <- Other.Mutations[-which(Other.Mutations == "C-18-T")]
  } else if (This.Mutation == "C-18-T") {
    Other.Mutations <- Other.Mutations[-which(Other.Mutations == "C-18-G")]
  }
  
  
  # open empty vectors that we'll fill in later with information
  P.Values.Vector <- vector()
  Magnitude.Effect.Vector <- vector()
  Mains.Vector <- vector()
  Epistatics.Vector <- vector()
  
  
  # loop through the other mutations with which This.Mutation could potentially interact
  for (i in 1:length(Other.Mutations)){
    
    print(i)
    
    # id of the mutation with which we'll perform a test to see if there is epistasis
    Each.Other.Mutation <- Other.Mutations[i]
    
    # which rows in This.Mutation.DF contain the other mutation?
    Rows.With.Other.Mutation <- which(sapply(X = as.character(This.Mutation.DF$Genotype.Final),
                                             FUN = function(x){
                                               Mutations.Here <- strsplit(x, ";")[[1]]
                                               Each.Other.Mutation %in% Mutations.Here
                                             }))
    
    
    # Calculate the effect [ln(A)] of This.Mutation in all the backgrounds of This.Mutation.DF
    All.Effects <- Logit.Transform(0.01*This.Mutation.DF$Final.PSI) - Logit.Transform(0.01*This.Mutation.DF$Starting.PSI)
    
    # Build a vector saying whether each genotype in This.Mutation.DF contains the potentially epistatic mutation
    Contains.Other.Mutation <- rep("No", nrow(This.Mutation.DF))
    Contains.Other.Mutation[Rows.With.Other.Mutation] <- "Yes"
    
    # A vector repeating the ID of the Other.Mutation the number of rows in This.Mutation.DF
    Other.Mutation.Vector <- rep(Each.Other.Mutation, nrow(This.Mutation.DF))
    
    # a data frame containing the three vectors we have just built
    Temporary.DF <- data.frame(Effect = All.Effects,
                               Contains.Other.Mutation = Contains.Other.Mutation,
                               Other.Mutation = Other.Mutation.Vector)
    
    # to perform a test, take all the effects in the presence of the Other.Mutation
    Effects.With.Other.Mutation <- All.Effects[Rows.With.Other.Mutation]
    Effects.With.Other.Mutation <- Effects.With.Other.Mutation[complete.cases(Effects.With.Other.Mutation)]
    
    # and in the absence of the Other.Mutation
    Effects.Without.Other.Mutation <- All.Effects[-Rows.With.Other.Mutation]
    Effects.Without.Other.Mutation <- Effects.Without.Other.Mutation[complete.cases(Effects.Without.Other.Mutation)]
    
    # if we have more than 3 observations in each group, perform the test
    if (length(Effects.With.Other.Mutation) >= 3 & length(Effects.Without.Other.Mutation) >=3) {
      T.Test.Object <- t.test(Effects.With.Other.Mutation, Effects.Without.Other.Mutation)
      P.Value <- T.Test.Object$p.value
      Magnitude.Effect <- T.Test.Object$estimate[1] - T.Test.Object$estimate[2]
    } else {
      P.Value <- NA
      Magnitude.Effect <- NA
    }
    
    # add the values we have calculated to the corresponding vectors:
    P.Values.Vector <- c(P.Values.Vector, P.Value)
    Magnitude.Effect.Vector <- c(Magnitude.Effect.Vector, Magnitude.Effect)
    Mains.Vector <- c(Mains.Vector, This.Mutation)
    Epistatics.Vector <- c(Epistatics.Vector, Each.Other.Mutation)
  }
  
  # a vector with the position of This.Mutation
  Position.Mutation <- rep(as.numeric(strsplit(This.Mutation, "-")[[1]][2]), length(P.Values.Vector))
  
  # a vector with the position of each of the Other.Mutations
  Position.Episatic.Mutations <- as.numeric(sapply(strsplit(Other.Mutations, "-"),
                                                   function(x){
                                                     x[2]
                                                   }))
  
  # distance between both mutations
  Distance <- abs(Position.Mutation - Position.Episatic.Mutations)
  
  # minus log10 of the p values from before...
  Minus.Log10.P.Values <- -1*log10(P.Values.Vector)
  
  # a data frame with information about the tests done with This.Mutation
  P.Value.Temp.DF <- data.frame(Main = Mains.Vector,
                                Epistatic = Epistatics.Vector,
                                Position.Main = Position.Mutation,
                                Position.Epistatic = Position.Episatic.Mutations,
                                Distance = Distance, 
                                Magnitude = Magnitude.Effect.Vector,
                                P.Value = P.Values.Vector,
                                Minus.Log10.P.Value = Minus.Log10.P.Values)
  
  # Combine the data frame above with the data frame containing information about all the mutations
  P.Values.Distance.DF <- rbind(P.Values.Distance.DF, P.Value.Temp.DF)
  
}

# fdr-adjust all the 130 p values
P.Values.Distance.DF$FDR <- p.adjust(p = P.Values.Distance.DF$P.Value, method = "fdr")

# minus log 10 the fdr scores
P.Values.Distance.DF$Minus.Log10.FDR <- -1*log10(P.Values.Distance.DF$FDR)

# convert some factors into characters because that's easier to work with
P.Values.Distance.DF$Main <- as.character(P.Values.Distance.DF$Main)
P.Values.Distance.DF$Epistatic <- as.character(P.Values.Distance.DF$Epistatic)























# Plot test results





# a matrix of all possible pairwise combinations except C18G and C18T
All.Possible.Interactions <- combn(x = names(Final.Vs.Starting.PSI.List), m = 2)
All.Possible.Interactions <- All.Possible.Interactions[,-47]

# set the FDR threshold to call something 'significant'
Significance.Level <- 0.05

# set a seed because I'll use jitter in the plots below and I want this to be reproducible
set.seed(2)

# start an empty plot
par(pty="s")
plot(NULL,
     xlim = c(0,35),
     ylim = c(0,80),
     xlab = "Distance",
     ylab = "-log10 FDR",
     las = 1)

# for each possible interaction, confirm whether both tests were significant at the FDR threshold defined above
for (i in 1:ncol(All.Possible.Interactions)) {
  
  # what are the two mutations involved in this potential interaction?
  Mut.1 <- All.Possible.Interactions[1,i]
  Mut.2 <- All.Possible.Interactions[2,i]
  
  # what are the two indices in P.Values.Distance.DF matching this interaction?
  First.Index <- which(P.Values.Distance.DF$Main == Mut.1 & P.Values.Distance.DF$Epistatic == Mut.2)
  Second.Index <- which(P.Values.Distance.DF$Main == Mut.2 & P.Values.Distance.DF$Epistatic == Mut.1)
  
  # is this interaction significant?
  Significant <- P.Values.Distance.DF$FDR[First.Index] < Significance.Level & P.Values.Distance.DF$FDR[Second.Index] < Significance.Level
  
  # jitter the bar (so we can see overlapping bars)
  Jittered.X <- jitter(x = P.Values.Distance.DF$Distance[First.Index], factor = 10)
  Colour <- "black"
  
  # if the interaction is significant, colour the bar orange
  if (Significant){
    Colour <- "orange"
  }
  
  # draw the bar
  segments(x0 = Jittered.X,
           y0 = P.Values.Distance.DF$Minus.Log10.FDR[First.Index],
           x1 = Jittered.X,
           y1 = P.Values.Distance.DF$Minus.Log10.FDR[Second.Index],
           lwd = 3,
           col = Colour)
}

# draw a dashed red line at the Significance level used as a threshold
abline(h = -1*log10(Significance.Level),
       lwd = 2, lty = 3, col = "red")







































# prepare for violin plot

# These are the columns (in  All.Possible.Interactions) that we found were significant
Significant.Interactions <- c(8, 14, 28, 32, 34, 54, 62)

# start two empty vectors, which we will fill with the distance between mutations involved
# in significant (Strong) or not significant (Weak) interactions
Distances.Weak.Interactions <- vector()
Distances.Strong.Interactions <- vector()

# loop through each of the interactions and:
# 1. look at the distance between both mutations
# 2. put this distance inside the Weak or Strong.Interactions vector
for (i in 1:ncol(All.Possible.Interactions))  {
  
  # get one index for this interaction
  Mut.1 <- All.Possible.Interactions[1,i]
  Mut.2 <- All.Possible.Interactions[2,i]
  First.Index <- which(P.Values.Distance.DF$Main == Mut.1 & P.Values.Distance.DF$Epistatic == Mut.2)
  
  # distance between mutations
  This.Distance <- P.Values.Distance.DF$Distance[First.Index]
  
  # send to distances vector
  if (i %in% Significant.Interactions){
    Distances.Strong.Interactions <- c(Distances.Strong.Interactions, This.Distance)
  } else {
    Distances.Weak.Interactions <- c(Distances.Weak.Interactions, This.Distance)
  }
  
}






# load violin plot library
library(vioplot)

# plot
set.seed(3)
plot(NULL,
     xlim = c(0.5,2.5),
     ylim = c(0,35),
     axes = F,
     xlab = "",
     ylab = "")
vioplot(Distances.Strong.Interactions,
        Distances.Weak.Interactions,
        col = "white",
        ylim = c(0,35),
        names = c("Strong", "Weak"),
        add = T, border = c("gray50"), rectCol = c("gray50"))
par(new=T)
Scatterplot.X <- c(rep(1,length(Distances.Strong.Interactions)), rep(2, length(Distances.Weak.Interactions)))
Scatterplot.Y <- c(Distances.Strong.Interactions, Distances.Weak.Interactions)
plot(x = jitter(Scatterplot.X, factor = 0.5),
     y = jitter(Scatterplot.Y, factor = 0.5),
     xlim = c(0.5,2.5),
     ylim = c(0,35),
     axes = F,
     xlab = "",
     ylab = "Distance",
     pch = 19,
     col = "black",
     cex = 0.8)
axis(side = 1,
     at = c(1,2),
     labels = c("Strong interaction", "No interaction"))
axis(side = 2,
     las = 1)

