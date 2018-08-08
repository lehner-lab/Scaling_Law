
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

Brain.PSI.Table <- PSI.Estimates[,Brain.Sample.IDs]
Skin.PSI.Table <- PSI.Estimates[,Skin.Sample.IDs]

# calculate mean PSI in brain
Brain.Mean.PSIs <- apply(X = Brain.PSI.Table,
                         MARGIN = 1,
                         FUN = mean,
                         na.rm=T)

# calculate mean PSI in skin
Skin.Mean.PSIs <- apply(X = Skin.PSI.Table,
                        MARGIN = 1,
                        FUN = mean,
                        na.rm=T)

# indices with bad data in each of the two tissues
Brain.To.Remove <- which(is.na(Brain.Mean.PSIs))
Skin.To.Remove <- which(is.na(Skin.Mean.PSIs))

# indices with bad data in either tissue
Rows.To.Remove <- union(Brain.To.Remove, Skin.To.Remove)

# remove datapoints that only have NA's
Brain.Mean.PSIs <- Brain.Mean.PSIs[-Rows.To.Remove]
Skin.Mean.PSIs <- Skin.Mean.PSIs[-Rows.To.Remove]

Skin.Minus.Brain <- Skin.Mean.PSIs - Brain.Mean.PSIs

# exons that are more included in skin
Exons.Up <- data.frame(Brain = Brain.Mean.PSIs[which(Skin.Minus.Brain > 0)],
                       Skin = Skin.Mean.PSIs[which(Skin.Minus.Brain > 0)])

# exons that are less included in skin
Exons.Down <- data.frame(Brain = Brain.Mean.PSIs[which(Skin.Minus.Brain < 0)],
                         Skin = Skin.Mean.PSIs[which(Skin.Minus.Brain < 0)])

# bin the data
Exons.Up$Group <- findInterval(x = Exons.Up$Brain,
                               vec = seq(0,1,0.1),
                               rightmost.closed = T)
Exons.Up$Group <- factor(Exons.Up$Group, levels = 1:10)

Exons.Down$Group <- findInterval(x = Exons.Down $Brain,
                                 vec = seq(0,1,0.1),
                                 rightmost.closed = T)
Exons.Down$Group <- factor(Exons.Down$Group, levels = 1:10)




##############
## 2. Plots ##
##############

# library
library(ggplot2)


###############################################################
# PLOT
#

ggplot(data = Exons.Up,
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
  coord_cartesian(ylim = c(0,100)) + 
  ylab(expression(Delta*PSI)) +
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

ggplot(data = Exons.Down,
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
  coord_cartesian(ylim = c(-100,0)) + 
  ylab(expression(Delta*PSI)) +
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
