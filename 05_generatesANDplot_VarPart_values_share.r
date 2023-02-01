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

library(Hmsc)
load("models_thin_1000_samples_250_chains_4.Rdata")

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
VPr_useintensity = data.frame(t(VPr$vals)*100,m$TrData$use_intensity)
write.csv2(VPr_useintensity, "VPr_useintensity.csv")

###Generates Vap plot (Fig. 1)
library(dplyr)

varpart<-VPr_useintensity

####Sum climatic and landsat variables####
climate<-varpart[,"Chelsa150_bio6"]+varpart[,"Chelsa150_bio15"]
Landsat<-varpart[,"Band3_median15"]+varpart[,"Band4_median15"]+varpart[,"Band7_median15"]
varpart2<-cbind(climate, Landsat, varpart[,-c(1:5)])
coln<-c("clim","Landsat", "poly_logCat", "log.HAND_50", "DistRivlog10", "hab_type2", "Random_transect", "nUses")
#put the variables name as first column
colnames(varpart2) <- coln
varpart2$NE<-100-rowSums(varpart2[1:7])
varpart2<-varpart2[with(varpart2, order(nUses, NE)),]
varpart3<-t(varpart2)
#add the variable names

#####BARPLOT####
                                    

par(mar = c(3, 10, 1, 1), lwd=0.1)
col=c("#FF0000", "#FFAA00", "#00AAFF",  "#0000FF" ,"#00FF00", "#00FFAA"  ,"cornsilk", "white")
barplot(varpart3[-c(8:9),],col = col, las=2, horiz=2, cex.names = 0.8)




