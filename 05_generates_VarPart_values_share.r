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
  
VPr_useintensity = data.frame(m$TrData$use_intensity,t(VPr$vals))
  
