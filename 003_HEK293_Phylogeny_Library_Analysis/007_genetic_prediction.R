################
## LOAD FILES ##
################

# load Whole.Dataset
load("002_whole_dataset.RData")



########################
## GENERAL PREP STUFF ##
########################

# use Whole.Dataset to get the genotypes with low variance
Low.Noise.Genotypes <- as.character(Whole.Dataset$Mutation.IDs)[which(Whole.Dataset$SD < 10)]


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



###############
## FUNCTIONS ##
###############

# 1) Function that returns Root Mean Squared Error
rmse <- function(error){
  sqrt(mean(error^2))
}


# 2) function to plot observed vs predicted
plotObservedVsPredictions <- function(Pred, Obs){
  par(pty = "s")
  
  # plot
  plot(NULL,
       xlim = c(-0.5,2),
       ylim = c(0,1),
       xlab = "",
       ylab = "",
       axes = F)
  abline(0,1,col = "gray80", lwd = 3)
  par(new=T)
  plot(Pred,
       Obs,
       xlim = c(-50,200),
       ylim = c(0,100),
       pch = 19,
       cex = 0.5,
       las = 1,
       xlab = "Predicted PSI",
       ylab = "Observed PSI")
}


# 3) residuals plot function
residualsPlot <- function(Pred, Residuals){
  par(pty = "s")
  
  # plot
  plot(NULL,
       ylim = c(-100,100),
       xlim = c(0,1),
       xlab = "",
       ylab = "",
       axes = F)
  abline(h=0,col = "gray80", lwd = 3)
  par(new=T)
  plot(Pred,
       Residuals,
       ylim = c(-100,100),
       pch = 19,
       cex = 0.5,
       las = 1,
       xlab = "Predicted PSI",
       ylab = "Observed - Predicted",
       xaxt = "n",
       yaxt = "n")
  
  # draw axes
  axis(side = 2, las = 1)
  axis(side = 1)
  
  # loess regression between X and Y
  Loess.X <- Pred
  Loess.Y <- Residuals
  Loess.Regression <- loess(Loess.Y ~ Loess.X)
  
  # x and y coordinates for loess line
  Line.X <- ceiling(min(Loess.X)):floor(max(Loess.X))
  Loess.Model <- predict(object = Loess.Regression,
                         newdata = data.frame(Loess.X = Line.X),
                         se = T)
  Line.Y <- Loess.Model$fit
  
  # loess confidence interval
  polygon(x = c(Line.X, rev(Line.X)),
          y = c(Loess.Model$fit - qt(0.975,Loess.Model$df)*Loess.Model$se,
                rev(Loess.Model$fit + qt(0.975,Loess.Model$df)*Loess.Model$se)),
          col = rgb(1,0,0,0.2),
          border = NA)
  
  # loess line
  lines(x = Line.X,
        y = Line.Y,
        col = "red",
        lwd = 3)
  
}











########################################################################
## MODEL 1: EFFECT OF MUTATIONS ON ANCESTRAL BG, LINEAR + INDEPENDENT ##
########################################################################

# Build a model based only on the effects of a mutation on the ancestral background

# start an empty matrix containing 12 columns, one for each mutation that could
# occur in an exon. These mutations will be variables in our model, where each
# variable will take either the value 1 (mutation present) or 0 (not present)
Variables.To.Train.Model <- matrix(data = 0,
                                   ncol = length(Names.In.Sequence.Order),
                                   nrow = 0)
colnames(Variables.To.Train.Model) <- Names.In.Sequence.Order



# this model will be trained on the ancestral sequence + the 12 single mutants
# empty vector...
Rows.With.Training.Genotypes <- vector()

# add ancestral exon index to 'Rows.With.Training.Genotypes'
Ancestor.Index <- which(as.character(Whole.Dataset$Mutation.IDs) == "")
Rows.With.Training.Genotypes <- c(Rows.With.Training.Genotypes, Ancestor.Index)


# add the 12 single mutant indices
for (i in 1:length(Names.In.Sequence.Order)) {
  Index <- which(as.character(Whole.Dataset$Mutation.IDs) == Names.In.Sequence.Order[i])
  Rows.With.Training.Genotypes <- c(Rows.With.Training.Genotypes, Index)
}
names(Rows.With.Training.Genotypes) <- c("Ancestor", as.character(Names.In.Sequence.Order))






# Collect the Genotype IDs from ancestor + 12 singles in a vector
Training.Genotypes.IDs <- as.character(Whole.Dataset$Mutation.IDs)[Rows.With.Training.Genotypes]

# now fill in the 'Variables.To.Train.Model' matrix with 13 rows of 1's and 0's,
# corresponding to the ancestor and 12 single mutants
invisible(sapply(Training.Genotypes.IDs,
                 function(x){
                   
                   # what mutations are there in this genotype?
                   singles.here <- strsplit(x, ";")[[1]]
                   
                   # build a one-row matrix with 12 columns full with 0's
                   Temporary.Matrix <- matrix(data = 0,
                                              ncol = length(Names.In.Sequence.Order),
                                              nrow = 1)
                   colnames(Temporary.Matrix) <- Names.In.Sequence.Order
                   
                   # for each mutation present in the genotype, change the 0 for a 1
                   for (each.single in singles.here){
                     Temporary.Matrix[1,each.single] <- 1
                   }
                   
                   # append the one-row matrix to the bottom of 'Variables.To.Train.Model'
                   Variables.To.Train.Model <<- rbind(Variables.To.Train.Model, Temporary.Matrix)
                   
                 }))




# We will now build a linear model using 'Variables.To.Train.Model' as the
# training set

# organise the training set in a data frame
Model.DF <- cbind(data.frame(Y = Whole.Dataset$Mean[Rows.With.Training.Genotypes]),
                  Variables.To.Train.Model)

# change the name of the variable names and remove all hyphens (easier to work with)
colnames(Model.DF)[2:13] <- sapply(as.character(colnames(Model.DF)[2:13]),
                                   function(x){
                                     paste(strsplit(x, "-")[[1]],
                                           sep = "",
                                           collapse = "")
                                   })

# build the model
Model <- lm(Y ~ (.), data = Model.DF)




# Now have to see if this model predicts the data













































# start an empty matrix containing 12 columns, one for each mutation that could
# occur in an exon. These mutations will be variables in our model, where each
# variable will take either the value 1 (mutation present) or 0 (not present)
Variables.To.Test.Model <- matrix(data = 0,
                                  ncol = length(Names.In.Sequence.Order),
                                  nrow = 0)
colnames(Variables.To.Test.Model) <- Names.In.Sequence.Order




# now fill in the 'Variables.To.Test.Model' matrix with lots of rows of 1's and 0's,
# corresponding to the ancestor and 12 single mutants
invisible(sapply(as.character(Whole.Dataset$Mutation.IDs),
                 function(x){
                   
                   # what mutations are there in this genotype?
                   singles.here <- strsplit(x, ";")[[1]]
                   
                   # build a one-row matrix with 12 columns full with 0's
                   Temporary.Matrix <- matrix(data = 0,
                                              ncol = length(Names.In.Sequence.Order),
                                              nrow = 1)
                   colnames(Temporary.Matrix) <- Names.In.Sequence.Order
                   
                   # for each mutation present in the genotype, change the 0 with a 1
                   for (each.single in singles.here){
                     Temporary.Matrix[1,each.single] <- 1
                   }
                   
                   # append the one-row matrix to the bottom of 'Variables.To.Train.Model'
                   Variables.To.Test.Model <<- rbind(Variables.To.Test.Model, Temporary.Matrix)
                   
                 }))




# change the name of the variable names and remove all hyphens
colnames(Variables.To.Test.Model) <- sapply(as.character(colnames(Variables.To.Test.Model)),
                                            function(x){
                                              paste(strsplit(x,"-")[[1]],
                                                    sep = "",
                                                    collapse = "")
                                            })









# remove stuff with PSI above 100 or sd > 10
Observations.To.Remove <- which(Whole.Dataset$Mean > 100 | Whole.Dataset$SD > 10)

# subset 'Variables.To.Test.Model' so it only includes stuff with PSI < 100
New.Data <- Variables.To.Test.Model[-Observations.To.Remove,]
New.Data <- as.data.frame(New.Data)
PSIs.For.New.Data <- Whole.Dataset$Mean[-Observations.To.Remove]






# use our model to make predictions
Predictions <- predict(Model, newdata = New.Data)

# plot
plotObservedVsPredictions(Pred = Predictions, Obs = PSIs.For.New.Data)

# plot
residualsPlot(Pred = Predictions, Residuals = PSIs.For.New.Data - Predictions)

# error in our predictions
Errors <- PSIs.For.New.Data - Predictions
rmse(Errors)



# set rule that predictions between 0 and 100
Predictions[which(Predictions < 0)] <- 0
Predictions[which(Predictions > 100)] <- 100



# plot
plotObservedVsPredictions(Pred = Predictions, Obs = PSIs.For.New.Data)
residualsPlot(Pred = Predictions, Residuals = PSIs.For.New.Data - Predictions)

# error in our predictions
Errors <- PSIs.For.New.Data - Predictions
rmse(Errors)









# now without restricting by SD

# remove stuff with PSI above 100 but any sd
Observations.To.Remove <- which(Whole.Dataset$Mean > 100)

# subset 'Variables.To.Test.Model' so it only includes stuff with PSI < 100
New.Data <- Variables.To.Test.Model[-Observations.To.Remove,]
New.Data <- as.data.frame(New.Data)
PSIs.For.New.Data <- Whole.Dataset$Mean[-Observations.To.Remove]






# use our model to make predictions
Predictions <- predict(Model, newdata = New.Data)


plotObservedVsPredictions(Pred = Predictions, Obs = PSIs.For.New.Data)

residualsPlot(Pred = Predictions, Residuals = PSIs.For.New.Data - Predictions)



# error in our predictions
Errors <- PSIs.For.New.Data - Predictions
rmse(Errors)





# set rule that predictions between 0 and 100
Predictions[which(Predictions < 0)] <- 0
Predictions[which(Predictions > 100)] <- 100



# plot
plotObservedVsPredictions(Pred = Predictions, Obs = PSIs.For.New.Data)
residualsPlot(Pred = Predictions, Residuals = PSIs.For.New.Data - Predictions)

# error in our predictions
Errors <- PSIs.For.New.Data - Predictions
rmse(Errors)



























##################################################################################
## 2. Model using average effect of mutations across many different backgrounds ##
##################################################################################





# New variables to train model -> use everything

# start an empty matrix containing 12 columns, one for each mutation that could
# occur in an exon. These mutations will be variables in our model, where each
# variable will take either the value 1 (mutation present) or 0 (not present)
Variables.To.Train.Model <- matrix(data = 0,
                                   ncol = length(Names.In.Sequence.Order),
                                   nrow = 0)
colnames(Variables.To.Train.Model) <- Names.In.Sequence.Order





# Collect the Genotype IDs from ancestor + 12 singles in a vector
Training.Genotypes.IDs <- as.character(Whole.Dataset$Mutation.IDs)

# now fill in the 'Variables.To.Train.Model' matrix with 13 rows of 1's and 0's,
# corresponding to the ancestor and 12 single mutants
invisible(sapply(Training.Genotypes.IDs,
                 function(x){
                   
                   # what mutations are there in this genotype?
                   singles.here <- strsplit(x, ";")[[1]]
                   
                   # build a one-row matrix with 12 columns full with 0's
                   Temporary.Matrix <- matrix(data = 0,
                                              ncol = length(Names.In.Sequence.Order),
                                              nrow = 1)
                   colnames(Temporary.Matrix) <- Names.In.Sequence.Order
                   
                   # for each mutation present in the genotype, change the 0 for a 1
                   for (each.single in singles.here){
                     Temporary.Matrix[1,each.single] <- 1
                   }
                   
                   # append the one-row matrix to the bottom of 'Variables.To.Train.Model'
                   Variables.To.Train.Model <<- rbind(Variables.To.Train.Model, Temporary.Matrix)
                   
                 }))



# We will now build a linear model using 'Variables.To.Train.Model' as the
# training set

# organise the training set in a data frame
Model.DF <- cbind(data.frame(Y = Whole.Dataset$Mean),
                  Variables.To.Train.Model)

# change the name of the variable names and remove all hyphens (easier to work with)
colnames(Model.DF)[2:13] <- sapply(as.character(colnames(Model.DF)[2:13]),
                                   function(x){
                                     paste(strsplit(x, "-")[[1]],
                                           sep = "",
                                           collapse = "")
                                   })







# remove stuff with PSI above 100 or sd > 10
Rows.To.Keep <- which(Whole.Dataset$SD < 10 & Whole.Dataset$Mean <= 100)

# build the model
Model <- lm(Y ~ (.), data = Model.DF[Rows.To.Keep,])












# subset 'Variables.To.Test.Model' so it only includes stuff with PSI < 100
New.Data <- Variables.To.Test.Model[Rows.To.Keep,]
New.Data <- as.data.frame(New.Data)
PSIs.For.New.Data <- Whole.Dataset$Mean[Rows.To.Keep]







# use our model to make predictions
Predictions <- predict(Model, newdata = New.Data)

#plot
plotObservedVsPredictions(Pred = Predictions, Obs = PSIs.For.New.Data)

residualsPlot(Pred = Predictions, Residuals = PSIs.For.New.Data - Predictions)


# error in our predictions


#########################################################
#  To 10-fold cross-validation & calculate RMSE #########

# 1. Set seed

# library with createFolds function
library(caret)

# function to calculate 10-fold cross validation
Ten.Fold.CV <- function(Model.Dataframe, Model.Formula, hardBounds = F, A.Space = F){
  
  # Create 10 folds
  Testing.Sets <- createFolds(Model.Dataframe$Y,
                              k = 10,
                              returnTrain = F)
  
  # a vector where we'll store the different RMSE values calcualted
  Vector.of.RMSEs <- c()
  
  # Calculate model RMSE for each of the 10 folds
  for (i in 1:10){
    
    # take a training and a testing set
    This.Training.Set <- Model.Dataframe[ -Testing.Sets[[i]], ]
    This.Testing.Set <- Model.Dataframe[ Testing.Sets[[i]], ]
    
    # Build a model for this fold
    Trained.Model <- lm(Model.Formula, data = This.Training.Set)
    
    # Make predictions
    Predictions.on.Testing.Set <- predict(object = Trained.Model, newdata = This.Testing.Set)
    
    # if there are hard bounds to be set, set them
    if (hardBounds) {
      if (length(which(Predictions.on.Testing.Set > 100)) > 0) {
        Predictions.on.Testing.Set[which(Predictions.on.Testing.Set > 100)] <- 100
      }
      
      if (length(which(Predictions.on.Testing.Set < 0)) > 0) {
        Predictions.on.Testing.Set[which(Predictions.on.Testing.Set < 0)] <- 0
      }
    }
    
    if (A.Space) {
      Predictions.on.Testing.Set <- 100*ConvertBackToPSI(Predictions.on.Testing.Set)
      # Calculate the error
      Errors <- 100*ConvertBackToPSI(This.Testing.Set$Y) - Predictions.on.Testing.Set
    } else {
      # Calculate the error
      Errors <- This.Testing.Set$Y - Predictions.on.Testing.Set
    }
    
    
    
    # calculate the rmse
    RMSE <- rmse(Errors)
    
    # add this rmse to the vector of rmses
    Vector.of.RMSEs <- c(Vector.of.RMSEs, RMSE)
    
  }
  
  # To combine RMSEs, use the formula in the best answer in:
  # https://stats.stackexchange.com/questions/85507/what-is-the-rmse-of-k-fold-cross-validation
  
  Overall.RMSE <- sqrt((sum(Vector.of.RMSEs^2))/10)
  print(Overall.RMSE)
  
}

# model formula
Model.Formula <- Y ~ .

# set the seed
set.seed(123)

# calculate 10-fold cross validation
Ten.Fold.CV(Model.Dataframe = Model.DF[Rows.To.Keep,],
            Model.Formula = Model.Formula)






# set rule that predictions between 0 and 100
Predictions[which(Predictions < 0)] <- 0
Predictions[which(Predictions > 100)] <- 100



# plot
plotObservedVsPredictions(Pred = Predictions, Obs = PSIs.For.New.Data)
residualsPlot(Pred = Predictions, Residuals = PSIs.For.New.Data - Predictions)

# set the seed
set.seed(123)

# calculate 10-fold cross validation
Ten.Fold.CV(Model.Dataframe = Model.DF[Rows.To.Keep,],
            Model.Formula = Model.Formula,
            hardBounds = T)




























###### without limiting by SD



# remove stuff with PSI above 100 or sd > 10
Rows.To.Keep <- which(Whole.Dataset$Mean <= 100)

# build the model
Model <- lm(Y ~ (.), data = Model.DF[Rows.To.Keep,])

# subset 'Variables.To.Test.Model' so it only includes stuff with PSI < 100
New.Data <- Variables.To.Test.Model[Rows.To.Keep,]
New.Data <- as.data.frame(New.Data)
PSIs.For.New.Data <- Whole.Dataset$Mean[Rows.To.Keep]


# use our model to make predictions
Predictions <- predict(Model, newdata = New.Data)

#plot
plotObservedVsPredictions(Pred = Predictions, Obs = PSIs.For.New.Data)

residualsPlot(Pred = Predictions, Residuals = PSIs.For.New.Data - Predictions)


# set the seed
set.seed(123)

# calculate 10-fold cross validation
Ten.Fold.CV(Model.Dataframe = Model.DF[Rows.To.Keep,],
            Model.Formula = Model.Formula)















# set rule that predictions between 0 and 100
Predictions[which(Predictions < 0)] <- 0
Predictions[which(Predictions > 100)] <- 100



# plot
plotObservedVsPredictions(Pred = Predictions, Obs = PSIs.For.New.Data)
residualsPlot(Pred = Predictions, Residuals = PSIs.For.New.Data - Predictions)

# set the seed
set.seed(123)

# calculate 10-fold cross validation
Ten.Fold.CV(Model.Dataframe = Model.DF[Rows.To.Keep,],
            Model.Formula = Model.Formula,
            hardBounds = T)




















































##########################################################################
## MODEL 3
##########################################################################


# remove stuff with PSI above 100 or sd > 10
Rows.To.Keep <- which(Whole.Dataset$SD < 10 & Whole.Dataset$Mean <= 100)

# build the model
Model <- lm(Y ~  . + C41G:C39T + T49C:G51C + T24C:G26T + T19G:C18G + T19G:C18T + C32T:G35T + C39T:G44A, data = Model.DF[Rows.To.Keep,])


# subset 'Variables.To.Test.Model' so it only includes stuff with PSI < 100
New.Data <- Variables.To.Test.Model[Rows.To.Keep,]
New.Data <- as.data.frame(New.Data)
PSIs.For.New.Data <- Whole.Dataset$Mean[Rows.To.Keep]


# use our model to make predictions
Predictions <- predict(Model, newdata = New.Data)

#plot
plotObservedVsPredictions(Pred = Predictions, Obs = PSIs.For.New.Data)

residualsPlot(Pred = Predictions, Residuals = PSIs.For.New.Data - Predictions)

# set the seed
set.seed(123)
Model.Formula <- Y ~  . + C41G:C39T + T49C:G51C + T24C:G26T + T19G:C18G + T19G:C18T + C32T:G35T + C39T:G44A
# calculate 10-fold cross validation
Ten.Fold.CV(Model.Dataframe = Model.DF[Rows.To.Keep,],
            Model.Formula = Model.Formula)




# set rule that predictions between 0 and 100
Predictions[which(Predictions < 0)] <- 0
Predictions[which(Predictions > 100)] <- 100



# plot
plotObservedVsPredictions(Pred = Predictions, Obs = PSIs.For.New.Data)
residualsPlot(Pred = Predictions, Residuals = PSIs.For.New.Data - Predictions)

# set the seed
set.seed(123)

# calculate 10-fold cross validation
Ten.Fold.CV(Model.Dataframe = Model.DF[Rows.To.Keep,],
            Model.Formula = Model.Formula,
            hardBounds = T)











######## without SD limitation
# remove stuff with PSI above 100
Rows.To.Keep <- which(Whole.Dataset$Mean <= 100)

# build the model
Model <- lm(Y ~  . + C41G:C39T + T49C:G51C + T24C:G26T + T19G:C18G + T19G:C18T + C32T:G35T + C39T:G44A, data = Model.DF[Rows.To.Keep,])


# subset 'Variables.To.Test.Model' so it only includes stuff with PSI < 100
New.Data <- Variables.To.Test.Model[Rows.To.Keep,]
New.Data <- as.data.frame(New.Data)
PSIs.For.New.Data <- Whole.Dataset$Mean[Rows.To.Keep]


# use our model to make predictions
Predictions <- predict(Model, newdata = New.Data)

#plot
plotObservedVsPredictions(Pred = Predictions, Obs = PSIs.For.New.Data)

residualsPlot(Pred = Predictions, Residuals = PSIs.For.New.Data - Predictions)

# set the seed
set.seed(123)
Model.Formula <- Y ~  . + C41G:C39T + T49C:G51C + T24C:G26T + T19G:C18G + T19G:C18T + C32T:G35T + C39T:G44A
# calculate 10-fold cross validation
Ten.Fold.CV(Model.Dataframe = Model.DF[Rows.To.Keep,],
            Model.Formula = Model.Formula)




# set rule that predictions between 0 and 100
Predictions[which(Predictions < 0)] <- 0
Predictions[which(Predictions > 100)] <- 100



# plot
plotObservedVsPredictions(Pred = Predictions, Obs = PSIs.For.New.Data)
residualsPlot(Pred = Predictions, Residuals = PSIs.For.New.Data - Predictions)

# set the seed
set.seed(123)

# calculate 10-fold cross validation
Ten.Fold.CV(Model.Dataframe = Model.DF[Rows.To.Keep,],
            Model.Formula = Model.Formula,
            hardBounds = T)






























































###############################################################
## MODEL 4
###############################################################



# function to calculate mutation effect A (PSIs in 0-1 scale)
CalculateA <- function(New.PSI, Old.PSI = 0.96) {
  log((New.PSI - Old.PSI*New.PSI) / (Old.PSI - Old.PSI*New.PSI))
}

# function to convert parameter A back into PSI (scale 0-1)
ConvertBackToPSI <- function(Parameter, Old.PSI = 0.96) {
  (exp(Parameter)*Old.PSI) / (1 - Old.PSI + exp(Parameter)*Old.PSI)
}




# organise the training set in a data frame
Model.DF <- cbind(data.frame(Y = sapply(X = Whole.Dataset$Mean/100,
                                        FUN = CalculateA)),
                  Variables.To.Train.Model)

# change the name of the variable names and remove all hyphens (easier to work with)
colnames(Model.DF)[2:13] <- sapply(as.character(colnames(Model.DF)[2:13]),
                                   function(x){
                                     paste(strsplit(x, "-")[[1]],
                                           sep = "",
                                           collapse = "")
                                   })

# remove stuff with PSI above 100 or sd > 10
Rows.To.Keep <- which(complete.cases(Model.DF) & Whole.Dataset$SD < 10)


# build the model
Model <- lm(Y ~ (.), data = Model.DF[Rows.To.Keep,])

















# subset 'Variables.To.Test.Model' so it only includes stuff with PSI < 100
New.Data <- Variables.To.Test.Model[Rows.To.Keep,]
New.Data <- as.data.frame(New.Data)
PSIs.For.New.Data <- Whole.Dataset$Mean[Rows.To.Keep]


# use our model to make predictions
Predictions <- predict(Model, newdata = New.Data)

#plot
plotObservedVsPredictions(Pred = 100*ConvertBackToPSI(Predictions), Obs = PSIs.For.New.Data)

residualsPlot(Pred = 100*ConvertBackToPSI(Predictions), Residuals = PSIs.For.New.Data - 100*ConvertBackToPSI(Predictions))









# set the seed
set.seed(123)
Model.Formula <- Y ~  . 
# calculate 10-fold cross validation
Ten.Fold.CV(Model.Dataframe = Model.DF[Rows.To.Keep,],
            Model.Formula = Model.Formula,
            A.Space = T)













######## without SD limitation
# remove stuff with PSI above 100
Rows.To.Keep <- which(Whole.Dataset$Mean <= 100)

# build the model
Model <- lm(Y ~  . , data = Model.DF[Rows.To.Keep,])

# subset 'Variables.To.Test.Model' so it only includes stuff with PSI < 100
New.Data <- Variables.To.Test.Model[Rows.To.Keep,]
New.Data <- as.data.frame(New.Data)
PSIs.For.New.Data <- Whole.Dataset$Mean[Rows.To.Keep]


# use our model to make predictions
Predictions <- predict(Model, newdata = New.Data)

#plot
plotObservedVsPredictions(Pred = 100*ConvertBackToPSI(Predictions),
                          Obs = PSIs.For.New.Data)
residualsPlot(Pred = 100*ConvertBackToPSI(Predictions),
              Residuals = PSIs.For.New.Data - 100*ConvertBackToPSI(Predictions))

# set the seed
set.seed(123)

# model formula
Model.Formula <- Y ~ . 

# calculate 10-fold cross validation rmse
Ten.Fold.CV(Model.Dataframe = Model.DF[Rows.To.Keep,],
            Model.Formula = Model.Formula,
            A.Space = T)









































































###########################################################################
#### MODEL 5
###########################################################################

# remove stuff with PSI above 100 or sd > 10
Rows.To.Keep <- which(Whole.Dataset$SD < 10 & Whole.Dataset$Mean <= 100)

# build the model
Model <- lm(Y ~  . + C41G:C39T + T49C:G51C + T24C:G26T + T19G:C18G + T19G:C18T + C32T:G35T + C39T:G44A, data = Model.DF[Rows.To.Keep,])


# subset 'Variables.To.Test.Model' so it only includes stuff with PSI < 100
New.Data <- Variables.To.Test.Model[Rows.To.Keep,]
New.Data <- as.data.frame(New.Data)
PSIs.For.New.Data <- Whole.Dataset$Mean[Rows.To.Keep]


# use our model to make predictions
Predictions <- predict(Model, newdata = New.Data)

#plot
plotObservedVsPredictions(Pred = 100*ConvertBackToPSI(Predictions), Obs = PSIs.For.New.Data)

residualsPlot(Pred = 100*ConvertBackToPSI(Predictions), Residuals = PSIs.For.New.Data - 100*ConvertBackToPSI(Predictions))

# set the seed
set.seed(123)
Model.Formula <- Y ~  . + C41G:C39T + T49C:G51C + T24C:G26T + T19G:C18G + T19G:C18T + C32T:G35T + C39T:G44A
# calculate 10-fold cross validation
Ten.Fold.CV(Model.Dataframe = Model.DF[Rows.To.Keep,],
            Model.Formula = Model.Formula, A.Space = T)















######## without SD limitation
# remove stuff with PSI above 100
Rows.To.Keep <- which(Whole.Dataset$Mean <= 100)

# build the model
Model <- lm(Y ~  . + C41G:C39T + T49C:G51C + T24C:G26T + T19G:C18G + T19G:C18T + C32T:G35T + C39T:G44A, data = Model.DF[Rows.To.Keep,])


# subset 'Variables.To.Test.Model' so it only includes stuff with PSI < 100
New.Data <- Variables.To.Test.Model[Rows.To.Keep,]
New.Data <- as.data.frame(New.Data)
PSIs.For.New.Data <- Whole.Dataset$Mean[Rows.To.Keep]


# use our model to make predictions
Predictions <- predict(Model, newdata = New.Data)

#plot
plotObservedVsPredictions(Pred = 100*ConvertBackToPSI(Predictions), Obs = PSIs.For.New.Data)

residualsPlot(Pred = 100*ConvertBackToPSI(Predictions), Residuals = PSIs.For.New.Data - 100*ConvertBackToPSI(Predictions))

# set the seed
set.seed(123)
Model.Formula <- Y ~  . + C41G:C39T + T49C:G51C + T24C:G26T + T19G:C18G + T19G:C18T + C32T:G35T + C39T:G44A
# calculate 10-fold cross validation
Ten.Fold.CV(Model.Dataframe = Model.DF[Rows.To.Keep,],
            Model.Formula = Model.Formula,
            A.Space = T)



