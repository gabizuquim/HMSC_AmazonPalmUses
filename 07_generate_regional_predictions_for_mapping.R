#################################################################################################################################################
#Codes associated to the analysis used the paper 
#"The legacy of human use in Amazonian palm communities along environmental and accessibility gradients". 
#Global Ecology and Biogeography. In press.
#Codes by Otso Ovaskainen, Mirkka Jones and Gabriela Zuquim
#The codes were modified from November 2020 Hmsc course scripts. 
#Current versions of these training materials are available at https://www.helsinki.fi/en/researchgroups/statistical-ecology/software/hmscCodes
#see readme file for detail
##################################################################################################################################################

#Import environmental covariate values for all grid cells across the Amazonian study region.
grid = fread("gridded_WAamazonian_predictorsV2.csv")

# Thin the Amazonian regional environmental predictor grid by taking every fourth x and y coordinate in the original grid 
# before making Hmsc predictions of species occurrence probabilities per cell computationally feasible.

uniqueX = sort(unique(grid$x))
uniqueY = sort(unique(grid$y))

thinX = seq(from = 1, to = length(uniqueX), by = 4)
selX = unique(grid$x)[thinX]
grid = grid[grid$x %in% selX,]

thinY = seq(from = 1, to = length(uniqueY), by = 4)
selY = unique(grid$y)[thinY]
grid = grid[grid$y %in% selY,]

fwrite(grid, "sparse_WAamazonian_predictors.csv")

##################

# NOTE: If no grid thinning is needed, predictions script starts here:

grid = fread("sparse_WAamazonian_predictors.csv")

library(Hmsc)
library(data.table)
library(abind)

#The fitted model (models_thin_1000_samples_250_chains_4.Rdata) 
#can be downloaded from https://etsin.fairdata.fi/dataset/46f59d1b-d925-4cf8-893e-37e50742a291

load("nonspatial_models_thin_1000_samples_250_chains_4.Rdata")
m = models[[1]]

m$XData = m$XData[, c("Chelsa150_bio6", "Chelsa150_bio15", "Band3_median15", "Band4_median15", "Band7_median15", "soilfern", "log.HAND_50")]

#Predictions are generated row by row (by Y coordinate) across the Amazonian regional predictor grid.
#WARNING! The next lines will take some hours to run in a regular computer
uniqueY = unique(grid$y)

predictions = data.frame(matrix(nrow = nrow(grid), ncol = ncol(m$Y)+2))
names(predictions) = c("x", "y", colnames(m$Y))

for(n in 1:length(uniqueY))
{
  selY = uniqueY[n]
  gridrow = grid[which(grid$y==selY),]
  
  xy.grid = data.frame(gridrow[, c("x", "y")])
  row.names(xy.grid) = gridrow$Pixel
  
  XData.grid = data.frame(gridrow[,-c(1:3)])
  
  #Check if some pixels cannot be predicted due to missing values:
  NAs = which(is.na(rowSums(XData.grid)))
  if(length(NAs>0))
  {
    xy.grid = xy.grid[-NAs,]
    XData.grid = XData.grid[-NAs,]
  }
  
  # Generate Hmsc predictions of probability of occurrence per species (predY = 1000 posterior predicted Y matrices, EpredY is their mean)
  Gradient = prepareGradient(m, XDataNew = XData.grid, sDataNew = list(transect = xy.grid))
  predY = predict(m, Gradient = Gradient, expected = TRUE)
  EpredY = apply(abind(predY,along = 3), c(1,2), mean)
  EpredY = cbind(xy.grid, EpredY)
  
  if(length(NAs>0))
  {
    NApixels = data.frame(cbind(gridrow$x[NAs], gridrow$y[NAs]))
    row.names(NApixels) = gridrow$Pixel[NAs]
    NApixels[,3:ncol(EpredY)] = NA
    names(NApixels) = names(EpredY)
    EpredY = rbind(NApixels, EpredY)
    EpredY = EpredY[order(as.numeric(row.names(EpredY))),]
  }
  
  predictions = rbind(predictions, EpredY)
  
}

# Based on the mean predicted occurrence probabilities per species, calculate predicted species richness (predS)
# and predicted trait means (predT) per grid cell.

preds = as.matrix(predictions[, 3:80])
predS = apply(predictions[, 3:80], 1, sum, na.rm = T)
predT = (preds %*% m$Tr)/matrix(rep(rowSums(preds),m$nt), ncol = m$nt)

#Combine posterior mean predicted species richness, human use trait values and species' p(occurrences):
predictions = cbind(predS, predT[,-1], predictions)

# Save data for mapping (Fig. 3), etc.
fwrite(predictions, "Palm_Predicted_Occurrence_Probabilities_WA.csv")
