#################################################################################################################################################
#Codes associated to the analysis used the paper 
#"The legacy of human use in Amazonian palm communities along environmental and accessibility gradients". 
#Global Ecology and Biogeography. In press.
#Codes by Otso Ovaskainen, Mirkka Jones and Gabriela Zuquim
#The codes were modified from November 2020 Hmsc course scripts. 
#Current versions of these training materials are available at https://www.helsinki.fi/en/researchgroups/statistical-ecology/software/hmscCodes
#see readme file for detail
##################################################################################################################################################
#The fitted model (models_thin_1000_samples_250_chains_4.Rdata) 
#can be downloaded from https://etsin.fairdata.fi/dataset/46f59d1b-d925-4cf8-893e-37e50742a291

load("models_thin_1000_samples_250_chains_4.Rdata")

library(Hmsc)

m = models[[1]]
postBeta = getPostEstimate(m, parName="Beta")
Betameans = t(postBeta$mean[-1,])
colnames(Betameans) = m$covNames[-1]
Fig2_data = cbind(m$TrData$use_intensity, Betameans)
colnames(Fig2_data)[1] = "N human uses"
colnames(Fig2_data)[-1] = paste0("Beta_", colnames(Fig2_data)[-1])

write.csv2(Fig2_data, "Fig2_data_UseIntensity_Beta.csv")
