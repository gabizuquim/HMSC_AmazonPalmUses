##add metadata####
library(Hmsc)

load(file = "unfitted_models.Rdata")

#### WARNING! The next lines can take several weeks or even months to run in a normal computer
#### For testing the codes purposes we suggest changing the values nchains, samples and thin to 2, 50, 1, respectively
#### Model fitting with the below settings took almost 2 weeks to run in a supercomputing center

nChains = 4
samples = 250
thin = 1000
for(n in 1:1)
{
  m = models[[n]]
  m = sampleMcmc(m, samples = samples, thin=thin,
                 adaptNf=rep(ceiling(0.4*samples*thin),m$nr),
                 transient = ceiling(0.5*samples*thin),
                 nChains = nChains, nParallel = 1)
  models[[n]] = m
}

filename_out = paste("models_thin_", as.character(thin),
                     "_samples_", as.character(samples),
                     "_chains_",as.character(nChains),
                     ".Rdata",sep = "")

save(models, modelnames, file=filename_out)
