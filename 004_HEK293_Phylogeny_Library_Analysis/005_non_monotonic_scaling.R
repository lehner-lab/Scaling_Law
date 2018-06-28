# load Whole.Dataset
load("002_whole_dataset.RData")

# use Whole.Dataset to get the genotypes with low variance
Low.Noise.Genotypes <- as.character(Whole.Dataset$Mutation.IDs)[which(Whole.Dataset$SD < 10)]

# load Final.Vs.Starting.PSI
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



###############################################################
# PLOT
#

par(mar = c(0,0,1,1) + 0.1)
par(oma = c(5,4,0,0) + 0.1)
par(mfrow=c(3,4))
par(pty="s")

for (i in 1:12){
  This.Number <- i
  
  # what is the ID of this mutation?
  This.Mutation <- as.character(Names.In.Sequence.Order)[This.Number]
  
  # extract the data frame
  This.Mutation.DF <- Final.Vs.Starting.PSI.List[[This.Mutation]]
  
  # Get low noise rows
  Low.Noise.Rows <- which(as.character(This.Mutation.DF$Genotype.Final) %in% Low.Noise.Genotypes & as.character(This.Mutation.DF$Genotype.Starting) %in% Low.Noise.Genotypes)
  
  # subset This.Mutation.DF
  This.Mutation.DF <- This.Mutation.DF[Low.Noise.Rows,]
  
  # plot final vs starting PSI
  par(pty="s")
  plot(NULL,
       xlim = c(0,100),
       ylim = c(0,100),
       xlab = "",
       ylab = "",
       axes = F)
  abline(0,1, col = "gray80", lwd = 2)
  par(new=T)
  plot(This.Mutation.DF$Starting.PSI,
       This.Mutation.DF$Final.PSI,
       pch = 19,
       cex = 0.5,
       xlim = c(0,100),
       ylim = c(0,100),
       xlab = "",
       ylab = "",
       main = "",
       axes = F)
  
  # surround the plot with a box (no axes were drawn yet)
  box()
  
  # stick the name of the mutation in the top left part of the plot
  text(x = 0,
       y = 90,
       labels = paste(strsplit(This.Mutation, split = "-")[[1]],
                      sep = "",
                      collapse = ""),
       pos = 4,
       cex = 1.5,
       col = "#3375CA")
  
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
title(xlab = "Starting PSI",
      ylab = "Final PSI",
      outer = TRUE,
      line = 3)

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
  
  # what is the ID of this mutation?
  This.Mutation <- as.character(Names.In.Sequence.Order)[This.Number]
  
  # extract the data frame
  This.Mutation.DF <- Final.Vs.Starting.PSI.List[[This.Mutation]]
  
  # plot final vs starting PSI
  par(pty="s")
  plot(NULL,
       xlim = c(0,100),
       ylim = c(0,100),
       xlab = "",
       ylab = "",
       axes = F)
  abline(0,1, col = "gray80", lwd = 2)
  par(new=T)
  plot(This.Mutation.DF$Starting.PSI,
       This.Mutation.DF$Final.PSI,
       pch = 19,
       cex = 0.5,
       xlim = c(0,100),
       ylim = c(0,100),
       xlab = "",
       ylab = "",
       main = "",
       axes = F)
  
  # surround the plot with a box (no axes were drawn yet)
  box()
  
  # stick the name of the mutation in the top left part of the plot
  text(x = 0,
       y = 90,
       labels = paste(strsplit(This.Mutation, split = "-")[[1]],
                      sep = "",
                      collapse = ""),
       pos = 4,
       cex = 1.5,
       col = "#3375CA")
  
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
title(xlab = "Starting PSI",
      ylab = "Final PSI",
      outer = TRUE,
      line = 3)

par(mar = c(5,4,4,2)+0.1)
par(oma = c(0,0,0,0))
par(mfrow=c(1,1))
        
#
#
###############################################################


# some mutations appear to display two distinct behaviours even before
# performing an epistasis analysis. If I'm going to fit curves to the
# data, it makes sense to fit different curves to each of these behaviours
# so make a note of what these interactions are:
Important.Epistatic.Pairs <- data.frame(Main =      c("C-18-G","T-19-G","C-32-T","T-49-C","G-51-C"),
                                        Epistatic = c("T-19-G","C-18-G","G-35-T","G-51-C","T-49-C"))


###############################################################
# PLOT
#

par(mar = c(0,0,1,1) + 0.1)
par(oma = c(5,4,0,0) + 0.1)
par(mfrow=c(3,4))
par(pty="s")

for (i in 1:12){
  This.Number <- i
  
  # what is the ID of this mutation?
  This.Mutation <- as.character(Names.In.Sequence.Order)[This.Number]
  
  # extract the data frame
  This.Mutation.DF <- Final.Vs.Starting.PSI.List[[This.Mutation]]
  
  # Get low noise rows
  Low.Noise.Rows <- which(as.character(This.Mutation.DF$Genotype.Final) %in% Low.Noise.Genotypes & as.character(This.Mutation.DF$Genotype.Starting) %in% Low.Noise.Genotypes)
  
  # subset This.Mutation.DF
  This.Mutation.DF <- This.Mutation.DF[Low.Noise.Rows,]
  
  # plot final vs starting PSI
  par(pty="s")
  plot(NULL,
       xlim = c(0,100),
       ylim = c(0,100),
       xlab = "",
       ylab = "",
       axes = F)
  abline(0,1, col = "gray80", lwd = 2)
  par(new=T)
  plot(This.Mutation.DF$Starting.PSI,
       This.Mutation.DF$Final.PSI,
       pch = 19,
       cex = 0.5,
       xlim = c(0,100),
       ylim = c(0,100),
       xlab = "",
       ylab = "",
       main = "",
       axes = F)
  
  # surround the plot with a box (no axes were drawn yet)
  box()
  
  # to fit the curves, we first need to check whether we'll fit one or two curves
  if (This.Mutation %in% as.character(Important.Epistatic.Pairs$Main)) {
    
    # if we fit 2 curves, check which is the epistatic mutation that causes the second behaviour
    Epistatic.Pair <- as.character(Important.Epistatic.Pairs$Epistatic)[which(as.character(Important.Epistatic.Pairs$Main) == This.Mutation)]
    
    # Which rows in This.Mutation.DF contain the epistatic mutation?
    Rows.With.Epistatic.Partner <- which(sapply(as.character(This.Mutation.DF$Genotype.Final),
                                                function(x){
                                                  Epistatic.Pair %in% strsplit(x, ";")[[1]]
                                                }))
    
    # Prepare to draw the first curve:
    # a) which starting + final genotypes contain the epistatic mutation?
    X <- 0.01*This.Mutation.DF$Starting.PSI[Rows.With.Epistatic.Partner]
    Y <- 0.01*This.Mutation.DF$Final.PSI[Rows.With.Epistatic.Partner]
    
    # b) logit transform
    Logit.X <- log((X)/(1-X))
    Logit.Y <- log((Y)/(1-Y))
    
    # c) stuff estimated > 100 will give an na, so need to remove those
    Remove.Indices <- union(which(is.na(Logit.X)), which(is.na(Logit.Y)))
    Logit.X <- Logit.X[-Remove.Indices]
    Logit.Y <- Logit.Y[-Remove.Indices]
    
    # d) calculate the effect ln(A) of our mutation
    Effect.K <- mean(Logit.Y - Logit.X)
    
    # e) a function for the relationship between Starting and Final PSI
    Curve.Fit <- function(X, K){
      (exp(K)*X)/(1-X+exp(K)*X)
    }
    
    # f) create the curve
    New.X <- seq(0,1,0.01)
    New.Y <- Curve.Fit(X = New.X, K = Effect.K)
    
    # g) and draw it
    lines(x = New.X*100,
          y = New.Y*100,
          col = "#3375CA",
          lwd = 2)
    
    # Now draw the second curve. Same as above, with the rows we didn't use earlier
    # a) which starting + final genotypes do NOT contain the epistatic mutation?
    X <- 0.01*This.Mutation.DF$Starting.PSI[-Rows.With.Epistatic.Partner]
    Y <- 0.01*This.Mutation.DF$Final.PSI[-Rows.With.Epistatic.Partner]
    
    # b) logit transform
    Logit.X <- log((X)/(1-X))
    Logit.Y <- log((Y)/(1-Y))
    
    # c) stuff estimated > 100 will give an na, so need to remove those
    Remove.Indices <- union(which(is.na(Logit.X)), which(is.na(Logit.Y)))
    Logit.X <- Logit.X[-Remove.Indices]
    Logit.Y <- Logit.Y[-Remove.Indices]
    
    # d) calculate the effect ln(A) of our mutation
    Effect.K <- mean(Logit.Y - Logit.X)
    
    # e) a function for the relationship between Starting and Final PSI
    Curve.Fit <- function(X, K){
      (exp(K)*X)/(1-X+exp(K)*X)
    }
    
    # f) create the curve
    New.X <- seq(0,1,0.01)
    New.Y <- Curve.Fit(X = New.X, K = Effect.K)
    
    # g) and draw it
    lines(x = 100*New.X,
          y = 100*New.Y,
          col = "#3375CA",
          lwd = 2)
    
  } else {
    
    # If there is only one clear behaviour, then we don't need to worry about epistasis
    # a) starting + final psis
    X <- 0.01*This.Mutation.DF$Starting.PSI
    Y <- 0.01*This.Mutation.DF$Final.PSI
    
    # b) logit transform
    Logit.X <- log((X)/(1-X))
    Logit.Y <- log((Y)/(1-Y))
    
    # c) stuff estimated > 100 will give an na, so need to remove those
    Remove.Indices <- union(which(is.na(Logit.X)), which(is.na(Logit.Y)))
    Logit.X <- Logit.X[-Remove.Indices]
    Logit.Y <- Logit.Y[-Remove.Indices]
    
    # d) calculate the effect ln(A) of our mutation
    Effect.K <- mean(Logit.Y - Logit.X)
    
    # e) a function for the relationship between Starting and Final PSI
    Curve.Fit <- function(X, K){
      (exp(K)*X)/(1-X+exp(K)*X)
    }
    
    # f) create the curve
    New.X <- seq(0,1,0.01)
    New.Y <- Curve.Fit(X = New.X, K = Effect.K)
    
    # g) and draw it
    lines(x = 100*New.X,
          y = 100*New.Y,
          col = "#3375CA",
          lwd = 2)
    
  }
  
  # stick the name of the mutation in the top left part of the plot
  text(x = 0,
       y = 90,
       labels = paste(strsplit(This.Mutation, split = "-")[[1]], sep = "", collapse = ""),
       pos = 4,
       cex = 1.5,
       col = "#00377F")
  
  # draw left-hand side axis if this plot is to the left
  if (i %in% c(1,5,9)){
    axis(side = 2, las = 1)
  }
  
  # draw bottom axis if this plot is at the bottom
  if (i %in% c(9,10,11,12)) {
    axis(side = 1)
  }
  
}

title(xlab = "Starting PSI",
      ylab = "Final PSI",
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
  
  # what is the ID of this mutation?
  This.Mutation <- as.character(Names.In.Sequence.Order)[This.Number]
  
  # extract the data frame
  This.Mutation.DF <- Final.Vs.Starting.PSI.List[[This.Mutation]]
  
  # Get low noise rows
  Low.Noise.Rows <- which(as.character(This.Mutation.DF$Genotype.Final) %in% Low.Noise.Genotypes & as.character(This.Mutation.DF$Genotype.Starting) %in% Low.Noise.Genotypes)
  
  # plot final vs starting PSI
  par(pty="s")
  plot(NULL,
       xlim = c(0,100),
       ylim = c(0,100),
       xlab = "",
       ylab = "",
       axes = F)
  abline(0,1, col = "gray80", lwd = 2)
  par(new=T)
  plot(This.Mutation.DF$Starting.PSI,
       This.Mutation.DF$Final.PSI,
       pch = 19,
       cex = 0.5,
       xlim = c(0,100),
       ylim = c(0,100),
       xlab = "",
       ylab = "",
       main = "",
       axes = F)
  
  # surround the plot with a box (no axes were drawn yet)
  box()
  
  # to fit the curves, we first need to check whether we'll fit one or two curves
  if (This.Mutation %in% as.character(Important.Epistatic.Pairs$Main)) {
    
    # if we fit 2 curves, check which is the epistatic mutation that causes the second behaviour
    Epistatic.Pair <- as.character(Important.Epistatic.Pairs$Epistatic)[which(as.character(Important.Epistatic.Pairs$Main) == This.Mutation)]
    
    # Which rows in This.Mutation.DF contain the epistatic mutation?
    Rows.With.Epistatic.Partner <- which(sapply(as.character(This.Mutation.DF$Genotype.Final)[Low.Noise.Rows],
                                                function(x){
                                                  Epistatic.Pair %in% strsplit(x, ";")[[1]]
                                                }))
    
    # Prepare to draw the first curve:
    # a) which starting + final genotypes contain the epistatic mutation?
    X <- 0.01*This.Mutation.DF$Starting.PSI[Low.Noise.Rows][Rows.With.Epistatic.Partner]
    Y <- 0.01*This.Mutation.DF$Final.PSI[Low.Noise.Rows][Rows.With.Epistatic.Partner]
    
    # b) logit transform
    Logit.X <- log((X)/(1-X))
    Logit.Y <- log((Y)/(1-Y))
    
    # c) stuff estimated > 100 will give an na, so need to remove those
    Remove.Indices <- union(which(is.na(Logit.X)), which(is.na(Logit.Y)))
    Logit.X <- Logit.X[-Remove.Indices]
    Logit.Y <- Logit.Y[-Remove.Indices]
    
    # d) calculate the effect ln(A) of our mutation
    Effect.K <- mean(Logit.Y - Logit.X)
    
    # e) a function for the relationship between Starting and Final PSI
    Curve.Fit <- function(X, K){
      (exp(K)*X)/(1-X+exp(K)*X)
    }
    
    # f) create the curve
    New.X <- seq(0,1,0.01)
    New.Y <- Curve.Fit(X = New.X, K = Effect.K)
    
    # g) and draw it
    lines(x = New.X*100,
          y = New.Y*100,
          col = "#3375CA",
          lwd = 2)
    
    # Now draw the second curve. Same as above, with the rows we didn't use earlier
    # a) which starting + final genotypes do NOT contain the epistatic mutation?
    X <- 0.01*This.Mutation.DF$Starting.PSI[Low.Noise.Rows][-Rows.With.Epistatic.Partner]
    Y <- 0.01*This.Mutation.DF$Final.PSI[Low.Noise.Rows][-Rows.With.Epistatic.Partner]
    
    # b) logit transform
    Logit.X <- log((X)/(1-X))
    Logit.Y <- log((Y)/(1-Y))
    
    # c) stuff estimated > 100 will give an na, so need to remove those
    Remove.Indices <- union(which(is.na(Logit.X)), which(is.na(Logit.Y)))
    Logit.X <- Logit.X[-Remove.Indices]
    Logit.Y <- Logit.Y[-Remove.Indices]
    
    # d) calculate the effect ln(A) of our mutation
    Effect.K <- mean(Logit.Y - Logit.X)
    
    # e) a function for the relationship between Starting and Final PSI
    Curve.Fit <- function(X, K){
      (exp(K)*X)/(1-X+exp(K)*X)
    }
    
    # f) create the curve
    New.X <- seq(0,1,0.01)
    New.Y <- Curve.Fit(X = New.X, K = Effect.K)
    
    # g) and draw it
    lines(x = 100*New.X,
          y = 100*New.Y,
          col = "#3375CA",
          lwd = 2)
    
  } else {
    
    # If there is only one clear behaviour, then we don't need to worry about epistasis
    # a) starting + final psis
    X <- 0.01*This.Mutation.DF$Starting.PSI[Low.Noise.Rows]
    Y <- 0.01*This.Mutation.DF$Final.PSI[Low.Noise.Rows]
    
    # b) logit transform
    Logit.X <- log((X)/(1-X))
    Logit.Y <- log((Y)/(1-Y))
    
    # c) stuff estimated > 100 will give an na, so need to remove those
    Remove.Indices <- union(which(is.na(Logit.X)), which(is.na(Logit.Y)))
    Logit.X <- Logit.X[-Remove.Indices]
    Logit.Y <- Logit.Y[-Remove.Indices]
    
    # d) calculate the effect ln(A) of our mutation
    Effect.K <- mean(Logit.Y - Logit.X)
    
    # e) a function for the relationship between Starting and Final PSI
    Curve.Fit <- function(X, K){
      (exp(K)*X)/(1-X+exp(K)*X)
    }
    
    # f) create the curve
    New.X <- seq(0,1,0.01)
    New.Y <- Curve.Fit(X = New.X, K = Effect.K)
    
    # g) and draw it
    lines(x = 100*New.X,
          y = 100*New.Y,
          col = "#3375CA",
          lwd = 2)
    
  }
  
  # stick the name of the mutation in the top left part of the plot
  text(x = 0,
       y = 90,
       labels = paste(strsplit(This.Mutation, split = "-")[[1]], sep = "", collapse = ""),
       pos = 4,
       cex = 1.5,
       col = "#00377F")
  
  # draw left-hand side axis if this plot is to the left
  if (i %in% c(1,5,9)){
    axis(side = 2, las = 1)
  }
  
  # draw bottom axis if this plot is at the bottom
  if (i %in% c(9,10,11,12)) {
    axis(side = 1)
  }
  
}

title(xlab = "Starting PSI",
      ylab = "Final PSI",
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
  
  # what is the ID of this mutation?
  This.Mutation <- as.character(Names.In.Sequence.Order)[This.Number]
  
  # extract the data frame
  This.Mutation.DF <- Final.Vs.Starting.PSI.List[[This.Mutation]]
  
  # Get low noise rows
  Low.Noise.Rows <- which(as.character(This.Mutation.DF$Genotype.Final) %in% Low.Noise.Genotypes & as.character(This.Mutation.DF$Genotype.Starting) %in% Low.Noise.Genotypes)
  
  # subset This.Mutation.DF
  This.Mutation.DF <- This.Mutation.DF[Low.Noise.Rows,]
  
  # plot final vs starting PSI
  par(pty="s")
  plot(NULL,
       xlim = c(0,100),
       ylim = c(-100,100),
       xlab = "",
       ylab = "",
       axes = F)
  abline(h=0, col = "gray80", lwd = 2)
  par(new=T)
  plot(This.Mutation.DF$Starting.PSI,
       This.Mutation.DF$Final.PSI - This.Mutation.DF$Starting.PSI,
       pch = 19,
       cex = 0.5,
       xlim = c(0,100),
       ylim = c(-100,100),
       xlab = "",
       ylab = "",
       main = "",
       axes = F)
  
  # surround the plot with a box (no axes were drawn yet)
  box()
  
  # to fit the curves, we first need to check whether we'll fit one or two curves
  if (This.Mutation %in% as.character(Important.Epistatic.Pairs$Main)) {
    
    # if we fit 2 curves, check which is the epistatic mutation that causes the second behaviour
    Epistatic.Pair <- as.character(Important.Epistatic.Pairs$Epistatic)[which(as.character(Important.Epistatic.Pairs$Main) == This.Mutation)]
    
    # Which rows in This.Mutation.DF contain the epistatic mutation?
    Rows.With.Epistatic.Partner <- which(sapply(as.character(This.Mutation.DF$Genotype.Final),
                                                function(x){
                                                  Epistatic.Pair %in% strsplit(x, ";")[[1]]
                                                }))
    
    # Prepare to draw the first curve:
    # a) which starting + final genotypes contain the epistatic mutation?
    X <- 0.01*This.Mutation.DF$Starting.PSI[Rows.With.Epistatic.Partner]
    Y <- 0.01*This.Mutation.DF$Final.PSI[Rows.With.Epistatic.Partner]
    
    # b) logit transform
    Logit.X <- log((X)/(1-X))
    Logit.Y <- log((Y)/(1-Y))
    
    # c) stuff estimated > 100 will give an na, so need to remove those
    Remove.Indices <- union(which(is.na(Logit.X)), which(is.na(Logit.Y)))
    Logit.X <- Logit.X[-Remove.Indices]
    Logit.Y <- Logit.Y[-Remove.Indices]
    
    # d) calculate the effect ln(A) of our mutation
    Effect.K <- mean(Logit.Y - Logit.X)
    
    # e) a function for the relationship between Starting and Final PSI
    Curve.Fit <- function(X, K){
      (exp(K)*X)/(1-X+exp(K)*X)
    }
    
    # f) create the curve
    New.X <- seq(0,1,0.01)
    New.Y <- Curve.Fit(X = New.X, K = Effect.K)
    
    # g) and draw it
    lines(x = New.X*100,
          y = New.Y*100 - New.X*100,
          col = "#3375CA",
          lwd = 2)
    
    # Now draw the second curve. Same as above, with the rows we didn't use earlier
    # a) which starting + final genotypes do NOT contain the epistatic mutation?
    X <- 0.01*This.Mutation.DF$Starting.PSI[-Rows.With.Epistatic.Partner]
    Y <- 0.01*This.Mutation.DF$Final.PSI[-Rows.With.Epistatic.Partner]
    
    # b) logit transform
    Logit.X <- log((X)/(1-X))
    Logit.Y <- log((Y)/(1-Y))
    
    # c) stuff estimated > 100 will give an na, so need to remove those
    Remove.Indices <- union(which(is.na(Logit.X)), which(is.na(Logit.Y)))
    Logit.X <- Logit.X[-Remove.Indices]
    Logit.Y <- Logit.Y[-Remove.Indices]
    
    # d) calculate the effect ln(A) of our mutation
    Effect.K <- mean(Logit.Y - Logit.X)
    
    # e) a function for the relationship between Starting and Final PSI
    Curve.Fit <- function(X, K){
      (exp(K)*X)/(1-X+exp(K)*X)
    }
    
    # f) create the curve
    New.X <- seq(0,1,0.01)
    New.Y <- Curve.Fit(X = New.X, K = Effect.K)
    
    # g) and draw it
    lines(x = 100*New.X,
          y = 100*New.Y - New.X*100,
          col = "#3375CA",
          lwd = 2)
    
  } else {
    
    # If there is only one clear behaviour, then we don't need to worry about epistasis
    # a) starting + final psis
    X <- 0.01*This.Mutation.DF$Starting.PSI
    Y <- 0.01*This.Mutation.DF$Final.PSI
    
    # b) logit transform
    Logit.X <- log((X)/(1-X))
    Logit.Y <- log((Y)/(1-Y))
    
    # c) stuff estimated > 100 will give an na, so need to remove those
    Remove.Indices <- union(which(is.na(Logit.X)), which(is.na(Logit.Y)))
    Logit.X <- Logit.X[-Remove.Indices]
    Logit.Y <- Logit.Y[-Remove.Indices]
    
    # d) calculate the effect ln(A) of our mutation
    Effect.K <- mean(Logit.Y - Logit.X)
    
    # e) a function for the relationship between Starting and Final PSI
    Curve.Fit <- function(X, K){
      (exp(K)*X)/(1-X+exp(K)*X)
    }
    
    # f) create the curve
    New.X <- seq(0,1,0.01)
    New.Y <- Curve.Fit(X = New.X, K = Effect.K)
    
    # g) and draw it
    lines(x = 100*New.X,
          y = 100*New.Y - New.X*100,
          col = "#3375CA",
          lwd = 2)
    
  }
  
  # stick the name of the mutation in the top left part of the plot
  text(x = 0,
       y = 80,
       labels = paste(strsplit(This.Mutation, split = "-")[[1]], sep = "", collapse = ""),
       pos = 4,
       cex = 1.5,
       col = "#00377F")
  
  # draw left-hand side axis if this plot is to the left
  if (i %in% c(1,5,9)){
    axis(side = 2, las = 1)
  }
  
  # draw bottom axis if this plot is at the bottom
  if (i %in% c(9,10,11,12)) {
    axis(side = 1)
  }
  
}

title(xlab = "Starting PSI",
      ylab = "Delta PSI",
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
  
  # what is the ID of this mutation?
  This.Mutation <- as.character(Names.In.Sequence.Order)[This.Number]
  
  # extract the data frame
  This.Mutation.DF <- Final.Vs.Starting.PSI.List[[This.Mutation]]
  
  # Get low noise rows
  Low.Noise.Rows <- which(as.character(This.Mutation.DF$Genotype.Final) %in% Low.Noise.Genotypes & as.character(This.Mutation.DF$Genotype.Starting) %in% Low.Noise.Genotypes)
  
  # plot final vs starting PSI
  par(pty="s")
  plot(NULL,
       xlim = c(0,100),
       ylim = c(-100,100),
       xlab = "",
       ylab = "",
       axes = F)
  abline(h=0, col = "gray80", lwd = 2)
  par(new=T)
  plot(This.Mutation.DF$Starting.PSI,
       This.Mutation.DF$Final.PSI - This.Mutation.DF$Starting.PSI,
       pch = 19,
       cex = 0.5,
       xlim = c(0,100),
       ylim = c(-100,100),
       xlab = "",
       ylab = "",
       main = "",
       axes = F)
  
  # surround the plot with a box (no axes were drawn yet)
  box()
  
  # to fit the curves, we first need to check whether we'll fit one or two curves
  if (This.Mutation %in% as.character(Important.Epistatic.Pairs$Main)) {
    
    # if we fit 2 curves, check which is the epistatic mutation that causes the second behaviour
    Epistatic.Pair <- as.character(Important.Epistatic.Pairs$Epistatic)[which(as.character(Important.Epistatic.Pairs$Main) == This.Mutation)]
    
    # Which rows in This.Mutation.DF contain the epistatic mutation?
    Rows.With.Epistatic.Partner <- which(sapply(as.character(This.Mutation.DF$Genotype.Final)[Low.Noise.Rows],
                                                function(x){
                                                  Epistatic.Pair %in% strsplit(x, ";")[[1]]
                                                }))
    
    # Prepare to draw the first curve:
    # a) which starting + final genotypes contain the epistatic mutation?
    X <- 0.01*This.Mutation.DF$Starting.PSI[Low.Noise.Rows][Rows.With.Epistatic.Partner]
    Y <- 0.01*This.Mutation.DF$Final.PSI[Low.Noise.Rows][Rows.With.Epistatic.Partner]
    
    # b) logit transform
    Logit.X <- log((X)/(1-X))
    Logit.Y <- log((Y)/(1-Y))
    
    # c) stuff estimated > 100 will give an na, so need to remove those
    Remove.Indices <- union(which(is.na(Logit.X)), which(is.na(Logit.Y)))
    Logit.X <- Logit.X[-Remove.Indices]
    Logit.Y <- Logit.Y[-Remove.Indices]
    
    # d) calculate the effect ln(A) of our mutation
    Effect.K <- mean(Logit.Y - Logit.X)
    
    # e) a function for the relationship between Starting and Final PSI
    Curve.Fit <- function(X, K){
      (exp(K)*X)/(1-X+exp(K)*X)
    }
    
    # f) create the curve
    New.X <- seq(0,1,0.01)
    New.Y <- Curve.Fit(X = New.X, K = Effect.K)
    
    # g) and draw it
    lines(x = New.X*100,
          y = New.Y*100 - New.X*100,
          col = "#3375CA",
          lwd = 2)
    
    # Now draw the second curve. Same as above, with the rows we didn't use earlier
    # a) which starting + final genotypes do NOT contain the epistatic mutation?
    X <- 0.01*This.Mutation.DF$Starting.PSI[Low.Noise.Rows][-Rows.With.Epistatic.Partner]
    Y <- 0.01*This.Mutation.DF$Final.PSI[Low.Noise.Rows][-Rows.With.Epistatic.Partner]
    
    # b) logit transform
    Logit.X <- log((X)/(1-X))
    Logit.Y <- log((Y)/(1-Y))
    
    # c) stuff estimated > 100 will give an na, so need to remove those
    Remove.Indices <- union(which(is.na(Logit.X)), which(is.na(Logit.Y)))
    Logit.X <- Logit.X[-Remove.Indices]
    Logit.Y <- Logit.Y[-Remove.Indices]
    
    # d) calculate the effect ln(A) of our mutation
    Effect.K <- mean(Logit.Y - Logit.X)
    
    # e) a function for the relationship between Starting and Final PSI
    Curve.Fit <- function(X, K){
      (exp(K)*X)/(1-X+exp(K)*X)
    }
    
    # f) create the curve
    New.X <- seq(0,1,0.01)
    New.Y <- Curve.Fit(X = New.X, K = Effect.K)
    
    # g) and draw it
    lines(x = 100*New.X,
          y = 100*New.Y - New.X*100,
          col = "#3375CA",
          lwd = 2)
    
  } else {
    
    # If there is only one clear behaviour, then we don't need to worry about epistasis
    # a) starting + final psis
    X <- 0.01*This.Mutation.DF$Starting.PSI[Low.Noise.Rows]
    Y <- 0.01*This.Mutation.DF$Final.PSI[Low.Noise.Rows]
    
    # b) logit transform
    Logit.X <- log((X)/(1-X))
    Logit.Y <- log((Y)/(1-Y))
    
    # c) stuff estimated > 100 will give an na, so need to remove those
    Remove.Indices <- union(which(is.na(Logit.X)), which(is.na(Logit.Y)))
    Logit.X <- Logit.X[-Remove.Indices]
    Logit.Y <- Logit.Y[-Remove.Indices]
    
    # d) calculate the effect ln(A) of our mutation
    Effect.K <- mean(Logit.Y - Logit.X)
    
    # e) a function for the relationship between Starting and Final PSI
    Curve.Fit <- function(X, K){
      (exp(K)*X)/(1-X+exp(K)*X)
    }
    
    # f) create the curve
    New.X <- seq(0,1,0.01)
    New.Y <- Curve.Fit(X = New.X, K = Effect.K)
    
    # g) and draw it
    lines(x = 100*New.X,
          y = 100*New.Y - New.X*100,
          col = "#3375CA",
          lwd = 2)
    
  }
    
  # stick the name of the mutation in the top left part of the plot
  text(x = 0,
       y = 80,
       labels = paste(strsplit(This.Mutation, split = "-")[[1]], sep = "", collapse = ""),
       pos = 4,
       cex = 1.5,
       col = "#00377F")
  
  # draw left-hand side axis if this plot is to the left
  if (i %in% c(1,5,9)){
    axis(side = 2, las = 1)
  }
  
  # draw bottom axis if this plot is at the bottom
  if (i %in% c(9,10,11,12)) {
    axis(side = 1)
  }
  
}

title(xlab = "Starting PSI",
      ylab = "Final PSI",
      outer = TRUE, line = 3)

par(mar = c(5,4,4,2)+0.1)
par(oma = c(0,0,0,0))
par(mfrow=c(1,1))
        
#
#
###############################################################
