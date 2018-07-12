
# library to read excel files
library(readxl)

# libraries for data wrangling
library(dplyr)
library(tidyr)

# read supplementary table 2
Ke.et.al.2018 <- read_excel("Data/All_HexMuts.xlsx")

# a unique genotype identifier
Ke.et.al.2018$ID <- rep(1:556,10)

# Build a table with the information we are interested
ES.Table <- Ke.et.al.2018[,c(1,2,8,9)]

# reshape to wide format
ES.Table <- ES.Table %>% spread(Hexmut, EI)
ES.Table <- as.data.frame(ES.Table)


###############################################################
# PLOT
#

# plot!
par(mfrow=c(3,3))
par(pty="s")
for (i in c(7,2,8,5,1,9,6,4,3)){
  plot(NULL,
       xlim = c(0,6),
       ylim = c(0,6),
       xlab = "",
       ylab = "",
       main = "",
       axes = F)
  abline(0,1, lwd = 2, col = "gray90")
  par(new=T)
  plot(ES.Table[which(ES.Table$Position > 23), 3],
       ES.Table[which(ES.Table$Position > 23), (3+i)],
       xlim = c(0,6),
       ylim = c(0,6),
       pch = 19,
       col = rgb(0,0,0,0.5),
       xlab = "Starting ES",
       ylab = "Final ES",
       las = 1,
       cex.lab = 1.5)
}

#
#
###############################################################
