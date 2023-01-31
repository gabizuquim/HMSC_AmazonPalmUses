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
VPr_useintensity = data.frame(m$TrData$use_intensity,t(VPr$vals)*100)
write.csv2(VPr_useintensity, "VPr_useintensity.csv")


###Generates Vap plot (Fig. 1)
library(dplyr)

varpart<-VPr_useintensity
####Sum climatic and landsat variables####
vpclim<-varpart[1:2,]
vpclim2<-colSums(vpclim[,-1])
vpRS<-varpart[3:5,]
vpRS2<-colSums(vpRS[,-1])

varpart2<-rbind(vpclim2, vpRS2,varpart[-c(1:5),-1])
#add the variable names
varpart2$X<-c("clim","Landsat", "poly_logCat", "log.HAND_50", "DistRivlog10", "hab_type2", "Random_transect", "TjurR2")
#put the variables name as first column
varpart2<-varpart2[,c(ncol(varpart2),1:(ncol(varpart2)-1))]
row.names(varpart2) <- varpart2$X 
varpart2 <- varpart2[,-1]
vp_prop<-t(varpart2)
vp_prop<-as.data.frame(vp_prop)
vp_prop$clim_p<-vp_prop$clim*vp_prop$TjurR2*100
vp_prop$Landsat_p<-vp_prop$Landsat*vp_prop$TjurR2*100
vp_prop$polylogCat_p<-vp_prop$poly_logCat*vp_prop$TjurR2*100
vp_prop$log.HAND_50_p<-vp_prop$log.HAND_50*vp_prop$TjurR2*100
vp_prop$DistRivlog10_p<-vp_prop$DistRivlog10*vp_prop$TjurR2*100
vp_prop$hab_type2_p<-vp_prop$hab_type2*vp_prop$TjurR2*100
vp_prop$Random_transect_p<-vp_prop$Random_transect*vp_prop$TjurR2*100
vp_prop$NE<-100-(vp_prop$TjurR2*100)
vp_prop<-t(vp_prop)
##just checking...
colSums(vp_prop[9:16,])
VP<-as.data.frame(vp_prop[c(11,9,10,12:16),])
VP<-t(VP)
VPuses=cbind(VP, use2)
VPuses2<-as.data.frame(VPuses[,c(1:8,16)])
VPuses2<-VPuses2[with(VPuses2, order(use_intensity, NE)),]
VPuses2<-t(VPuses2)
VPuses3<-VPuses2[1:8,]



#####BARPLOT####
                                    # Add legend to barplot
       legend = c("Soil (9.1)","Climate (11.7)","Landsat (4.3)","HAND (0.5)", "Distance to river (0.4)", 
                  "Habitat type (3.8)", "Random: transect (4.5)","Unexplained (0.6)")
       

pdf("barplot_uses.pdf", width = 10, height = 12)
par(mar = c(3, 10, 1, 1), lwd=0.1)
col=c("#FF0000", "#FFAA00", "#00AAFF",  "#0000FF" ,"#00FF00", "#00FFAA"  ,"cornsilk", "white")

barplot(VPuses3,col = col, las=2, horiz=2, cex.names = 0.8)
legend("bottomright",                                    # Add legend to barplot
       legend = legend,
       fill = col,
       bg="white")

dev.off()

pdf("LEGEND.pdf", width = 10, height = 12)
par(mar = c(3, 10, 1, 1), lwd=0.1)
col=c("#FF0000", "#FFAA00", "#00AAFF",  "#0000FF" ,"#00FF00", "#00FFAA"  ,"cornsilk", "white")

barplot(VPuses3,col = col, las=2, horiz=2, cex.names = 0.8)
legend("bottomright",                                    # Add legend to barplot
       legend = legend,
       fill = col,
       bg="white", cex=2)

dev.off()

#####ADD lines to create suppl Stable 1#####


