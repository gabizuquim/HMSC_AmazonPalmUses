#S4
#setwd("") # set directory to the folder where the folders "data", "models" and "panels" are
#setwd("/Volumes/HD_Elements/workspaceHD_Gabi_Zuquim25OCT2017/Danish project 2019/Ms.3 Palms HMSC/analysis/scripts_FULL_HMSCframework")
library(Hmsc)
install.packages("colorspace")
library(colorspace)
install.packages("vioplot")
library("vioplot")
install.packages("sm")
library(sm)

#include in samples_list and thin_list only those models that you have actually fitted!
#samples_list = c(5,250,250,250)
#thin_list = c(1,1,10,100)
#nst = length(thin_list)
#nChains = 4

#samples_list = c(5,20)
#thin_list = c(1,1)
#nChains = 2
nst = length(thin_list)

ma = NULL
na = NULL
for (Lst in 1:nst) {
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  
  
  filename = paste("models/models_thin_", as.character(thin),
                   "_samples_", as.character(samples),
                   "_chains_",as.character(nChains),".Rdata",sep = "")
  load(filename)
  nm = length(models)
  for(j in 1:nm){
    mpost = convertToCodaObject(models[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
    tmp = summary(psrf.beta)
    if(is.null(ma)){
      ma=psrf.beta[,1]
      na = paste0(as.character(thin),",",as.character(samples))
    } else {
      ma = cbind(ma,psrf.beta[,1])
      if(j==1){
        na = c(na,paste0(as.character(thin),",",as.character(samples)))
      } else {
        na = c(na,"")
      }
    }
  }
}
##the lines below are not working because I failed to install vioplot
pdf(file=paste("panels/MCMC_convergence.pdf"))
par(mfrow=c(2,1))
vioplot(ma,col=rainbow_hcl(nm),names=na,ylim=c(0,max(ma)),main="psrf(beta)")
vioplot(ma,col=rainbow_hcl(nm),names=na,ylim=c(0.9,1.1),main="psrf(beta)")
dev.off()

####Try vioplot with ggplot2
##the lines below are not working because I failed to install vioplot
library(ggplot2)
#pdf(file=paste("panels/MCMC_convergence_ggplot2.pdf"))
#par(mfrow=c(2,1))
#p<-ggplot(ma)+
#  geom_violin(aes(fill=rainbow_hcl(nm)))
              #,names=na,ylim=c(0,max(ma)),main="psrf(beta)"))
#vioplot(ma,col=rainbow_hcl(nm),names=na,ylim=c(0.9,1.1),main="psrf(beta)")
#dev.off()

