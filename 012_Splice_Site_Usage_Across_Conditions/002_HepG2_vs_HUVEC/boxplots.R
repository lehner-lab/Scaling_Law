
#################################
## Alternative 5' Splice Sites ##
#################################


Samples.Table <- read.table(file = "Data/Huvec_HepG2_TABLE_ALTD.txt",
                            sep = "\t",
                            header = TRUE)

# Calculate the mean PSU for all alternative splice site events in HepG2
Samples.Table$Mean.HepG2 <- apply(X = Samples.Table[,grep("^HepG2",
                                                          colnames(Samples.Table),
                                                          perl = T)],
                                  MARGIN = 1,
                                  FUN = function(x){
                                    psi.values <- as.numeric(x[c(TRUE, FALSE)])
                                    psi.qualities <- as.character(x[c(FALSE, TRUE)])
                                    
                                    values.to.use <- which(psi.qualities == "Pass")
                                    
                                    if (length(values.to.use) > 0){
                                      mean.psi <- mean(psi.values[values.to.use])
                                    } else {
                                      mean.psi <- NA
                                    }
                                    
                                    mean.psi
                                  })

# Calculate the mean PSU for all events in Huvec
Samples.Table$Mean.Huvec <- apply(X = Samples.Table[,grep("^Huvec",
                                                          colnames(Samples.Table),
                                                          perl = T)],
                                  MARGIN = 1,
                                  FUN = function(x){
                                    
                                    psi.values <- as.numeric(x[c(TRUE, FALSE)])
                                    psi.qualities <- as.character(x[c(FALSE, TRUE)])
                                    
                                    values.to.use <- which(psi.qualities == "Pass")
                                    
                                    if (length(values.to.use) > 0){
                                      mean.psi <- mean(psi.values[values.to.use])
                                    } else {
                                      mean.psi <- NA
                                    }
                                    
                                    mean.psi
                                  })

Samples.Table$Huvec.Minus.HepG2 <- Samples.Table$Mean.Huvec - Samples.Table$Mean.HepG2


PSU.Down <- Samples.Table[which(Samples.Table$Huvec.Minus.HepG2 < 0),] 
PSU.Up <- Samples.Table[which(Samples.Table$Huvec.Minus.HepG2 > 0),] 


PSU.Up$Group <- findInterval(x = PSU.Up$Mean.HepG2,
                             vec = seq(0,100,10),
                             rightmost.closed = T)
PSU.Up$Group <- factor(PSU.Up $Group,
                       levels = 1:10)

PSU.Down$Group <- findInterval(x = PSU.Down$Mean.HepG2,
                               vec = seq(0,100,10),
                               rightmost.closed = T)
PSU.Down$Group <- factor(PSU.Down$Group,
                         levels = 1:10)


# library for plotting stuff
library(ggplot2)

###############################################################
# PLOT
#

ggplot(data = PSU.Up, mapping = aes(x = Group,
                                    y = Huvec.Minus.HepG2)) +
  geom_boxplot(outlier.shape = NA,
               notch = T,
               fill = "#D66F79") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        aspect.ratio = 1,
        axis.text.x = element_text(size = 10,
                                   angle = 90,
                                   vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  coord_cartesian(ylim = c(0,50)) + 
  ylab(expression(Delta*PSU)) +
  xlab("Starting PSI") +
  scale_x_discrete(labels = c("1" = "[0-10)",
                              "2" = "[10-20)",
                              "3" = "[20-30)",
                              "4" = "[30-40)",
                              "5" = "[40-50)",
                              "6" = "[50-60)",
                              "7" = "[60-70)",
                              "8" = "[70-80)",
                              "9" = "[80-90)",
                              "10" = "[90-100]"))

#
#
###############################################################


###############################################################
# PLOT
#

ggplot(data = PSU.Down,
       mapping = aes(x = Group,
                     y = Huvec.Minus.HepG2)) +
  geom_boxplot(outlier.shape = NA,
               notch = T,
               fill = "#6EA7D3") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        aspect.ratio = 1,
        axis.text.x = element_text(size = 10,
                                   angle = 90,
                                   vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  coord_cartesian(ylim = c(-50,0)) + 
  ylab(expression(Delta*PSU)) +
  xlab("Starting PSI") +
  scale_x_discrete(labels = c("1" = "[0-10)",
                              "2" = "[10-20)",
                              "3" = "[20-30)",
                              "4" = "[30-40)",
                              "5" = "[40-50)",
                              "6" = "[50-60)",
                              "7" = "[60-70)",
                              "8" = "[70-80)",
                              "9" = "[80-90)",
                              "10" = "[90-100]"))
#
#
###############################################################










#################################
## Alternative 3' Splice Sites ##
#################################


Samples.Table <- read.table(file = "Data/Huvec_HepG2_TABLE_ALTA.txt",
                            sep = "\t",
                            header = TRUE)

# Calculate the mean PSU for all alternative splice site events in HepG2
Samples.Table$Mean.HepG2 <- apply(X = Samples.Table[,grep("^HepG2",
                                                          colnames(Samples.Table),
                                                          perl = T)],
                                  MARGIN = 1,
                                  FUN = function(x){
                                    psi.values <- as.numeric(x[c(TRUE, FALSE)])
                                    psi.qualities <- as.character(x[c(FALSE, TRUE)])
                                    
                                    values.to.use <- which(psi.qualities == "Pass")
                                    
                                    if (length(values.to.use) > 0){
                                      mean.psi <- mean(psi.values[values.to.use])
                                    } else {
                                      mean.psi <- NA
                                    }
                                    
                                    mean.psi
                                  })

# Calculate the mean PSU for all events in Huvec
Samples.Table$Mean.Huvec <- apply(X = Samples.Table[,grep("^Huvec",
                                                          colnames(Samples.Table),
                                                          perl = T)],
                                  MARGIN = 1,
                                  FUN = function(x){
                                    
                                    psi.values <- as.numeric(x[c(TRUE, FALSE)])
                                    psi.qualities <- as.character(x[c(FALSE, TRUE)])
                                    
                                    values.to.use <- which(psi.qualities == "Pass")
                                    
                                    if (length(values.to.use) > 0){
                                      mean.psi <- mean(psi.values[values.to.use])
                                    } else {
                                      mean.psi <- NA
                                    }
                                    
                                    mean.psi
                                  })

Samples.Table$Huvec.Minus.HepG2 <- Samples.Table$Mean.Huvec - Samples.Table$Mean.HepG2


PSU.Down <- Samples.Table[which(Samples.Table$Huvec.Minus.HepG2 < 0),] 
PSU.Up <- Samples.Table[which(Samples.Table$Huvec.Minus.HepG2 > 0),] 


PSU.Up$Group <- findInterval(x = PSU.Up$Mean.HepG2,
                             vec = seq(0,100,10),
                             rightmost.closed = T)
PSU.Up$Group <- factor(PSU.Up $Group,
                       levels = 1:10)

PSU.Down$Group <- findInterval(x = PSU.Down$Mean.HepG2,
                               vec = seq(0,100,10),
                               rightmost.closed = T)
PSU.Down$Group <- factor(PSU.Down$Group,
                         levels = 1:10)


# library for plotting stuff
library(ggplot2)

###############################################################
# PLOT
#

ggplot(data = PSU.Up, mapping = aes(x = Group,
                                    y = Huvec.Minus.HepG2)) +
  geom_boxplot(outlier.shape = NA,
               notch = T,
               fill = "#D66F79") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        aspect.ratio = 1,
        axis.text.x = element_text(size = 10,
                                   angle = 90,
                                   vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  coord_cartesian(ylim = c(0,50)) + 
  ylab(expression(Delta*PSU)) +
  xlab("Starting PSI") +
  scale_x_discrete(labels = c("1" = "[0-10)",
                              "2" = "[10-20)",
                              "3" = "[20-30)",
                              "4" = "[30-40)",
                              "5" = "[40-50)",
                              "6" = "[50-60)",
                              "7" = "[60-70)",
                              "8" = "[70-80)",
                              "9" = "[80-90)",
                              "10" = "[90-100]"))

#
#
###############################################################


###############################################################
# PLOT
#

ggplot(data = PSU.Down,
       mapping = aes(x = Group,
                     y = Huvec.Minus.HepG2)) +
  geom_boxplot(outlier.shape = NA,
               notch = T,
               fill = "#6EA7D3") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        aspect.ratio = 1,
        axis.text.x = element_text(size = 10,
                                   angle = 90,
                                   vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)) +
  coord_cartesian(ylim = c(-50,0)) + 
  ylab(expression(Delta*PSU)) +
  xlab("Starting PSI") +
  scale_x_discrete(labels = c("1" = "[0-10)",
                              "2" = "[10-20)",
                              "3" = "[20-30)",
                              "4" = "[30-40)",
                              "5" = "[40-50)",
                              "6" = "[50-60)",
                              "7" = "[60-70)",
                              "8" = "[70-80)",
                              "9" = "[80-90)",
                              "10" = "[90-100]"))
#
#
###############################################################