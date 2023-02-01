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
library(ggplot2)
library(abind)

#The fitted model (models_thin_1000_samples_250_chains_4.Rdata) 
#can be downloaded from https://etsin.fairdata.fi/dataset/46f59d1b-d925-4cf8-893e-37e50742a291

load("models_thin_1000_samples_250_chains_4.Rdata")
m = models[[1]]
covariates = all.vars(m$XFormula)[c(6, 8, 7, 1:5)]
gradientnames = c("log(Soil cations)", "log(Distance to river)", "HAND", "Min. temp. coldest month", "Precipitation seasonality", "Landsat Band 3", "Landsat Band 4", "Landsat Band 7")

col_richness = rgb(0.93,0.3,0.17)
col_food = "cyan"
col_material = rgb(0.67, 1, 0)
col_medicine = rgb(1, 0.8, 0)
col_intensity = "blue"

legendfill = c(col_food, col_medicine, col_material)

cicol_richness = adjustcolor(col_richness, alpha.f=0.1)
cicol_food = adjustcolor(col_food, alpha.f=0.1)
cicol_material = adjustcolor(col_material, alpha.f=0.1)
cicol_medicine = adjustcolor(col_medicine, alpha.f=0.15)
cicol_intensity = adjustcolor(col_intensity, alpha.f=0.1)

maxmin = function(x){max(x)>x[1]}

# Traits, Species richness vs gradients: Supplementary figure.

correlations = list()

pdf("SuppFig1_TraitsRichness_per_gradient.pdf", width = 5.5, height = 11)         
par(mfrow = c(8, 3), mar = c(4, 4, 1, 1), mgp = c(2, 0.75, 0))

for(n in 1:8)
{
  covariate = covariates[n]
  gradientname = gradientnames[n]
  Gradient2 = constructGradient(m,focalVariable = covariate,non.focalVariables = 1)
  focalgradient = Gradient2$XDataNew[,1]
  
  predY2 = predict(m, Gradient = Gradient2, expected = TRUE)
  
  predT = lapply(predY2, function(a) (a %*% m$Tr)/matrix(rep(rowSums(a),m$nt), ncol = m$nt))
  predT = abind(predT, along = 3)
  qpredT = apply(predT, MARGIN = c(1, 2), quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  predS = abind(lapply(predY2, rowSums),along=2)
  qpredS = apply(predS, c(1), quantile, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
  
  correlations[[n]] = round(cor(data.frame(qpredS[2,], qpredT[2,,])[, c(1, 3:6)]), 2)

  lo = qpredS[1, ]
  hi = qpredS[3, ]
  me = qpredS[2, ]
  
  linetype = 1
  linewidth = 1.5
                 
  plot(focalgradient, me, ylim  = c(0,28), xlab = gradientname, ylab = "Species richness", cex.lab = 0.85, cex.axis = 0.7, type = "l", col = col_richness, lwd = linewidth)
  polygon(c(focalgradient, rev(focalgradient)), c(lo, rev(hi)), col =  cicol_richness, border = FALSE)
  
  for(i in 1:4)
  {
    focaltrait = switch(i, "use_intensity", "human_food", "material", "medicine")
    cicol = switch(i, cicol_intensity, cicol_food, cicol_material, cicol_medicine)
    colour = switch(i, col_intensity, col_food, col_material, col_medicine)

            qpredT_focal = qpredT[, , focaltrait]
            lo = qpredT_focal[1, ]
            hi = qpredT_focal[3, ]
            me = qpredT_focal[2, ]
            
            #Pr = mean(apply(predT[, focaltrait, ], 2, maxmin))
            #linetype = ifelse(Pr>0.95|Pr<0.05, 1, 3)
          
            if(i == 1)
            {
              plot(focalgradient, me, ylim  = c(2,6), xlab = gradientname, ylab = "Use intensity", cex.lab = 0.85, cex.axis = 0.7, type = "l", lty = linetype, col = colour, lwd = linewidth)
              polygon(c(focalgradient, rev(focalgradient)), c(lo, rev(hi)), col =  cicol, border = FALSE)
            }
            if(i == 2)
            {
            plot(focalgradient, me, ylim  = c(0.4,1), xlab = gradientname, ylab = "Use proportion", cex.lab = 0.85, cex.axis = 0.7, type = "l", lty = linetype, col = colour, lwd = linewidth)
            polygon(c(focalgradient, rev(focalgradient)), c(lo, rev(hi)), col =  cicol, border = FALSE)
            }
            if(i > 2)
            {
              points(focalgradient, me, type = "l", col = colour, lwd = linewidth, lty = linetype)
              polygon(c(focalgradient, rev(focalgradient)), c(lo, rev(hi)), col =  cicol, border = FALSE)
            }
            if(i == 4)
            {
              legend("bottomleft", legend = c("F", "Me.", "Ma."), horiz = TRUE, bty = "n", border = FALSE, fill = legendfill, cex=0.7)
            }
  }
  }

dev.off()

rich_use_cortable = data.frame(correlations[[1]][-1,1], correlations[[2]][-1,1], correlations[[3]][-1,1], correlations[[4]][-1,1], correlations[[5]][-1,1], correlations[[6]][-1,1], correlations[[7]][-1,1], correlations[[8]][-1,1])
names(rich_use_cortable) = gradientnames

write.csv2(rich_use_cortable, "rich_use_cortable.csv", row.names = F)
