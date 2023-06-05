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
library(dplyr)

load("models_thin_1000_samples_250_chains_4.Rdata")

m = models[[1]]
preds = computePredictedValues(m)
VP = computeVariancePartitioning(m)
MF = evaluateModelFit(hM=m, predY=preds)
R2 = MF$TjurR2
for(k in 1:m$ns)
{VP$vals[,k] = R2[k]*VP$vals[,k]}


### Prepare data and plot Variation partitioning results (Fig. 1)

varpart = data.frame(t(VP$vals))
# Sum variation explained by climatic and landsat variables, respectively:
vpclim = rowSums(varpart[,1:2])
vplandsat = rowSums(varpart[,3:5])
varpart = cbind(varpart[,6], vpclim, vplandsat, varpart[,-c(1:6)], MF$TjurR2)
colnames(varpart) = c("poly_logCat", "Climate", "Landsat", "log.HAND_50", "DistRivlog10", "hab_type2", "Random_transect", "TjurR2")
varpart$Unexplained = 1-MF$TjurR2

#Convert values into percentages rather than proportions:
varpart = round(varpart*100,2)

#Add human use trait data per palm species:
VPuses = cbind(varpart, m$TrData[, c("human_food", "material", "medicine", "use_intensity")])
VPuses = VPuses[order(VPuses$use_intensity, VPuses$Unexplained),]

#Save the combined data as a csv file:
write.csv2(VPuses, "VP_TjurR2_Human_uses.csv")

#####BARPLOT####

#Subset VP values per species for plotting:
VPuses_to_plot<-t(VPuses)[c(1:7),]

#Check rowMeans for legend:
round(rowMeans(VPuses_to_plot),1)
round(mean(VPuses$Unexplained),1)

# Create legend for barplot
legend = c("Soil (9.1)","Climate (11.7)","Landsat (4.3)","HAND (0.5)", "Distance to river (0.4)", 
           "Habitat type (3.8)", "Random: transect (4.5)","Unexplained (65.6)")

pdf("Barplot_VP_useintensity.pdf", width = 10, height = 12)
par(mar = c(4, 10, 1, 1), lwd=0.1)
colPlot=c("#FF0000", "#FFAA00", "#00AAFF",  "#0000FF" ,"#00FF00", "#00FFAA", "cornsilk")
colLegend=c("#FF0000", "#FFAA00", "#00AAFF",  "#0000FF" ,"#00FF00", "#00FFAA"  ,"cornsilk", "white")

barplot(VPuses_to_plot,col = colPlot, las=2, horiz=2, cex.names = 0.8, xlim = c(0,100), xlab = "Variance proportion (%)")
legend("bottomright",
       legend = legend,
       fill = colLegend,
       cex=0.9)
textYcoord = c(4.1, 24, 35, 44, 49, 64.6, 77, 84)
for(n in 1:8)
{
  text(70, textYcoord[n], n-1, cex = 1.5) 
}

dev.off()




