#################################################################################################################################################
#Codes associated to the analysis used the paper 
#"The legacy of human use in Amazonian palm communities along environmental and accessibility gradients". 
#Global Ecology and Biogeography. In press.
#Codes by Otso Ovaskainen, Mirkka Jones and Gabriela Zuquim
#The codes were modified from November 2020 Hmsc course scripts. 
#Current versions of these training materials are available at https://www.helsinki.fi/en/researchgroups/statistical-ecology/software/hmscCodes
#see readme file for detail
##################################################################################################################################################

library(Hmsc)
#For computing-time optimization, cross-validation was run on a model fitted with 10 x fewer samples than the version used in the rest of the paper.
#To fit this model change the value of samples to 100 in script S3.
#Models with a large number of iterations are heavy to run. The output of the spatial model 
##fitted for the purpose of the paper associated with these code has been made open is a data repository.

load(file = "models_thin_100_samples_250_chains_4.Rdata") #object "models_thin_1000_samples_250_chains_4.Rdata" 
                                                          #can be accessed here: https://doi.org/10.23729/9e2da5a5-5c6c-45a6-b23a-e7af08aface2

mpost = convertToCodaObject(models[[1]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
    tmp = summary(psrf.beta)
    ma=psrf.beta[,1]
      
vioplot(ma,ylim=c(min(ma),max(ma)),main="psrf(beta)")



###WARNING!!! The below lines can take weeks to months to run in a regular computer
#for shorter runs, reduce nChains, samples, thin and nfolds to smaller values, e.g. 2, 50, 1, 2, respectively.
nChains = 4
samples = 250
thin = 100
nfolds = 4

MF = list()
MFCV = list()
WAIC = list()

  m = models[[1]]
  preds = computePredictedValues(m)
  MF[[1]] = evaluateModelFit(hM=m, predY=preds)
  partition = createPartition(m, nfolds = nfolds)
  preds = computePredictedValues(m, partition=partition)
  MFCV[[1]] = evaluateModelFit(hM=m, predY=preds)
  WAIC[[1]] = computeWAIC(m)       

filename_out = paste("MF_models_thin_", as.character(thin),
                     "_samples_", as.character(samples),
                     "_chains_",as.character(nChains),
                     ".Rdata",sep = "")

save(MF, MFCV, WAIC, modelnames, file = filename_out)
