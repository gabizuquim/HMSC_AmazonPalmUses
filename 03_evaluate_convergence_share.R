####METADATA####

library(Hmsc)
install.packages("colorspace")
library(colorspace)
install.packages("vioplot")
library("vioplot")
install.packages("sm")
library(sm)

load("models_thin_1000_samples_250_chains_4.Rdata")

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

par(mfrow=c(2,1))
vioplot(ma,col=rainbow_hcl(nm),names=na,ylim=c(0,max(ma)),main="psrf(beta)")
vioplot(ma,col=rainbow_hcl(nm),names=na,ylim=c(0.9,1.1),main="psrf(beta)")
