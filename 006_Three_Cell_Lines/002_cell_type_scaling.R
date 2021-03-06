####################
## 1. Hek to Hela ##
####################

# RT-PCR results for HEK293
Hek <- c(98.3,
         99.1,
         97.9,
         99,
         95.5,
         98.2,
         73.5,
         94.3,
         27.2,
         94.3,
         3.5,
         15.1,
         3.4,
         44)
Hek.SD <- c(0.2,
            0.7,
            0.8,
            0.2,
            0.3,
            0.2,
            0.7,
            0.5,
            1.6,
            0.5,
            1.5,
            0.9,
            1,
            1.5)

# RT-PCR results for HeLa
Hela <- c(98,
          98.4,
          96.3,
          95.3,
          95.2,
          97.4,
          89.6,
          97.6,
          59.3,
          97.9,
          8.7,
          55.1,
          11.1,
          82.6)
Hela.SD <- c(0.5,
             0.5,
             0.4,
             1.6,
             0.7,
             0.4,
             6.3,
             0.4,
             0.9,
             0.4,
             2.4,
             1.4,
             1.9,
             0.8)

# propagated error for delta PSI
Combined.SD <- sqrt(Hela.SD^2 + Hek.SD^2)


###############################################################
# PLOT
#

par(pty="s")
plot(NULL,
     xlim = c(0,100),
     ylim = c(0,100),
     xlab = "",
     ylab = "",
     axes = F)
abline(0,1, lwd = 3, col = "gray90")
segments(x0 = Hek,
         x1 = Hek,
         y0 = Hela+Hela.SD,
         y1 = Hela-Hela.SD,
         col = "gray70",
         lwd = 1.5)
segments(x0 = Hek+Hek.SD,
         x1 = Hek-Hek.SD,
         y0 = Hela,
         y1 = Hela,
         col = "gray70",
         lwd = 1.5)
par(new=T)
plot(Hek,
     Hela,
     xlim = c(0,100),
     ylim = c(0,100),
     pch = 19,
     las = 1,
     xlab = "HEK293",
     ylab = "HeLa")


# library for orthogonal least squares
library(onls)

# find the parameter A that best explains the data
Scaling.Effect <- onls(Hela ~ (100*A*Hek)/(100-Hek+A*Hek), start = list(A=1) )

# and draw the curve
Line.X <- seq(0,100,1)
Line.Y <- predict(object = Scaling.Effect,
                  newdata = data.frame(Hek = Line.X))
lines(x = Line.X,
      y = Line.Y,
      col = "orange",
      lwd = 2)
        
#
#
###############################################################


###############################################################
# PLOT
#

par(pty="s")
plot(NULL,
     xlim = c(0,100),
     ylim = c(-100,100),
     xlab = "",
     ylab = "",
     axes = F)
abline(h=0, lwd = 3, col = "gray90")
segments(x0 = Hek,
         x1 = Hek,
         y0 = Hela-Hek+Combined.SD,
         y1 = Hela-Hek-Combined.SD,
         col = "gray70",
         lwd = 1.5)
segments(x0 = Hek+Hek.SD,
         x1 = Hek-Hek.SD,
         y0 = Hela-Hek,
         y1 = Hela-Hek,
         col = "gray70",
         lwd = 1.5)
par(new=T)
plot(Hek,Hela-Hek,
     xlim = c(0,100),
     ylim = c(-100,100),
     pch = 19,
     las = 1,
     xlab = "HEK293",
     ylab = "HeLa - HEK293")

lines(x = Line.X,
      y = Line.Y-Line.X,
      col = "orange",
      lwd = 2)
        
#
#
###############################################################




######################
## 2. Hela to COS-7 ##
######################

# RT-PCR results for COS-7
Cos7 <- c(98.5,
          97.5,
          95.3,
          96.9,
          98.2,
          100,
          84.3,
          98.6,
          32.3,
          98.2,
          1.6,
          17.8,
          0.2,
          48)
Cos7.SD <- c(1.3,
             1.4,
             0.4,
             1,
             1.1,
             0.1,
             2.4,
             0.6,
             1.5,
             1.6,
             0.2,
             0.7,
             0.6,
             1.5)

# propagated error for delta PSI
Combined.SD <- sqrt(Hela.SD^2 + Cos7.SD^2)


###############################################################
# PLOT
#

par(pty="s")
plot(NULL,
     xlim = c(0,100),
     ylim = c(0,100),
     xlab = "",
     ylab = "",
     axes = F)
abline(0,1, lwd = 3, col = "gray90")
segments(x0 = Hela,
         x1 = Hela,
         y0 = Cos7+Cos7.SD,
         y1 = Cos7-Cos7.SD,
         col = "gray70",
         lwd = 1.5)
segments(x0 = Hela+Hela.SD,
         x1 = Hela-Hela.SD,
         y0 = Cos7,
         y1 = Cos7,
         col = "gray70",
         lwd = 1.5)
par(new=T)
plot(Hela,
     Cos7,
     xlim = c(0,100),
     ylim = c(0,100),
     pch = 19,
     las = 1,
     xlab = "HeLa",
     ylab = "COS-7")

# find the parameter A that best explains the data
Scaling.Effect <- onls(Cos7 ~ (100*A*Hela)/(100-Hela+A*Hela), start = list(A=1) )

# and draw the curve
Line.X <- seq(0,100,1)
Line.Y <- predict(object = Scaling.Effect,
                  newdata = data.frame(Hela = Line.X))
lines(x = Line.X,
      y = Line.Y,
      col = "orange",
      lwd = 2)
        
#
#
###############################################################


###############################################################
# PLOT
#

par(pty="s")
plot(NULL,
     xlim = c(0,100),
     ylim = c(-100,100),
     xlab = "",
     ylab = "",
     axes = F)
abline(h=0, lwd = 3, col = "gray90")
segments(x0 = Hela,
         x1 = Hela,
         y0 = Cos7-Hela+Combined.SD,
         y1 = Cos7-Hela-Combined.SD,
         col = "gray70",
         lwd = 1.5)
segments(x0 = Hela+Hela.SD,
         x1 = Hela-Hela.SD,
         y0 = Cos7-Hela,
         y1 = Cos7-Hela,
         col = "gray70",
         lwd = 1.5)
par(new=T)
plot(Hela,Cos7-Hela,
     xlim = c(0,100),
     ylim = c(-100,100),
     pch = 19,
     las = 1,
     xlab = "HeLa",
     ylab = "COS-7 - HeLa")

lines(x = Line.X,
      y = Line.Y-Line.X,
      col = "orange",
      lwd = 2)
        
#
#
###############################################################
