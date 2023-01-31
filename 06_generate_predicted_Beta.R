load("models_thin_1000_samples_250_chains_4.Rdata")

library(Hmsc)

m = models[[1]]
postBeta = getPostEstimate(m, parName="Beta")
Betameans = t(postBeta$mean[-1,])
colnames(Betameans) = m$covNames[-1]
head(m$TrData)
Fig2_data = cbind(m$TrData$use_intensity, Betameans)
colnames(Fig2_data)[1] = "N human uses"
colnames(Fig2_data)[-1] = paste0("Beta_", colnames(Fig2_data)[-1])

Write.csv2(Fig2_data, "Fig2_data_UseIntensity_Beta.csv")
