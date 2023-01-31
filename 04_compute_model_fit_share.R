###METADATA####
library(Hmsc)
#For computing-time optimization, cross-validation was run on a model fitted with 10 x fewer samples than the version used in the rest of the paper.
#To fit this model change the value of samples to 100 in script S3.
load(file = "models_thin_100_samples_250_chains_4.Rdata")

nChains = 4
samples = 250
thin = 100

MF = list()
MFCV = list()
WAIC = list()

for(n in 1:1){
  m = models[[n]]
  preds = computePredictedValues(m)
  MF[[n]] = evaluateModelFit(hM=m, predY=preds)
  partition = createPartition(m, nfolds = 4)
  preds = computePredictedValues(m, partition=partition, nParallel = nChains)
  MFCV[[n]] = evaluateModelFit(hM=m, predY=preds)
  WAIC[[n]] = computeWAIC(m)       
}

filename_out = paste("MF_models_thin_", as.character(thin),
                     "_samples_", as.character(samples),
                     "_chains_",as.character(nChains),
                     ".Rdata",sep = "")

save(MF, MFCV, WAIC, modelnames, file = filename_out)
