#METADATA

library(Hmsc)
load("models_thin_1000_samples_250_chains_4.Rdata")

library(Hmsc)

  m = models[[1]]
  preds = computePredictedValues(m)
  VP = computeVariancePartitioning(m)
  VPr = VP
  MF = evaluateModelFit(hM=m, predY=preds)
  R2 = MF$TjurR2
  for(k in 1:m$ns)
  {VPr$vals[,k] = R2[k]*VPr$vals[,k]}

#combine the variation partitioning values with use diversity trait data for each species for figure 1
#var part values are scaled by Tjur´s values
VPr_useintensity = data.frame(m$TrData$use_intensity,t(VPr$vals))
write.csv2(VPr_useintensity, "VPr_useintensity.csv")

#####ADD lines to create suppl Stable 1#####


