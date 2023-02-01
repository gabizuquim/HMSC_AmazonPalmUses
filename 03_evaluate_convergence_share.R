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
install.packages("colorspace")
library(colorspace)
install.packages("vioplot")
library("vioplot")
install.packages("sm")
library(sm)

load("models_thin_1000_samples_250_chains_4.Rdata")

mpost = convertToCodaObject(models[[1]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
    tmp = summary(psrf.beta)
    ma=psrf.beta[,1]
      
vioplot(ma,ylim=c(min(ma),max(ma)),main="psrf(beta)")

