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

#Models with a large number of iterations are heavy to run. The output of the spatial model 
##fitted for the purpose of the paper associated with these code has been made open is a data repository.

load("models_thin_1000_samples_250_chains_4.Rdata")#object "models_thin_1000_samples_250_chains_4.Rdata" can be accessed here: https://doi.org/10.23729/9e2da5a5-5c6c-45a6-b23a-e7af08aface2

mpost = convertToCodaObject(models[[1]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
    tmp = summary(psrf.beta)
    ma=psrf.beta[,1]
      
vioplot(ma,ylim=c(min(ma),max(ma)),main="psrf(beta)")

