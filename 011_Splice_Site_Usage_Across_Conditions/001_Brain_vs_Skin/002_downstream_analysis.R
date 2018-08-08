
#############################
## 1. Data processing in R ##
#############################

# load file
load("Tissues_Compared_Datasets.RData")

# how many samples for each tissue?
table(GTEX.Dataset$Tissue)

# which row numbers refer to a sample in brain or skin?
Brain.Indices <- which(GTEX.Dataset$Tissue == "Brain")
Skin.Indices <- which(GTEX.Dataset$Tissue == "Skin")

# use those row numbers as indices to extract relevant IDs
Brain.Sample.IDs <- as.character(GTEX.Dataset$ID)[Brain.Indices]
Skin.Sample.IDs <- as.character(GTEX.Dataset$ID)[Skin.Indices]

Brain.PSU.Table <- PSU.Estimates[,Brain.Sample.IDs]
Skin.PSU.Table <- PSU.Estimates[,Skin.Sample.IDs]

# calculate mean PSU in brain
Brain.Mean.PSUs <- apply(X = Brain.PSU.Table,
                         MARGIN = 1,
                         FUN = mean,
                         na.rm=T)

# calculate mean PSU in skin
Skin.Mean.PSUs <- apply(X = Skin.PSU.Table,
                        MARGIN = 1,
                        FUN = mean,
                        na.rm=T)

# indices with bad data in each of the two tissues
Brain.To.Remove <- which(is.na(Brain.Mean.PSUs))
Skin.To.Remove <- which(is.na(Skin.Mean.PSUs))

# indices with bad data in either tissue
Rows.To.Remove <- union(Brain.To.Remove, Skin.To.Remove)

# remove datapoints that only have NA's
Brain.Mean.PSUs <- Brain.Mean.PSUs[-Rows.To.Remove]
Skin.Mean.PSUs <- Skin.Mean.PSUs[-Rows.To.Remove]

Skin.Minus.Brain <- Skin.Mean.PSUs - Brain.Mean.PSUs

# exons that are more included in skin
PSU.Up <- data.frame(Brain = Brain.Mean.PSUs[which(Skin.Minus.Brain > 0)],
                     Skin = Skin.Mean.PSUs[which(Skin.Minus.Brain > 0)])

# exons that are less included in skin
PSU.Down <- data.frame(Brain = Brain.Mean.PSUs[which(Skin.Minus.Brain < 0)],
                       Skin = Skin.Mean.PSUs[which(Skin.Minus.Brain < 0)])

# bin the data
PSU.Up$Group <- findInterval(x = PSU.Up$Brain,
                             vec = seq(0,1,0.1),
                             rightmost.closed = T)
PSU.Up$Group <- factor(PSU.Up$Group, levels = 1:10)

PSU.Down$Group <- findInterval(x = PSU.Down$Brain,
                               vec = seq(0,1,0.1),
                               rightmost.closed = T)
PSU.Down$Group <- factor(Exons.Down$Group, levels = 1:10)




##############
## 2. Plots ##
##############

# library
library(ggplot2)


###############################################################
# PLOT
#

ggplot(data = PSU.Up,
       mapping = aes(x = Group,
                     y = (Skin - Brain)*100)) +
  geom_boxplot(outlier.shape = NA,
               notch = T,
               fill = "#D66F79") +
  # geom_jitter(aes(Group,(Skin - Brain)*100),
  #             position=position_jitter(width=0.25,
  #                                      height=0),
  #             alpha=0.1,
  #             size=0.5,
  #             show.legend=FALSE) +
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
  xlab("Starting PSU") +
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
                     y = (Skin - Brain)*100)) +
  geom_boxplot(outlier.shape = NA,
               notch = T,
               fill = "#6EA7D3") +
  # geom_jitter(aes(Group,(Skin - Brain)*100),
  #             position=position_jitter(width=0.25,
  #                                      height=0),
  #             alpha=0.1,
  #             size=0.5,
  #             show.legend=FALSE) +
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
  xlab("Starting PSU") +
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
