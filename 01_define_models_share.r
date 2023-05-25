#################################################################################################################################################
# Codes associated with the analyses presented in the paper 
# "The legacy of human use in Amazonian palm communities along environmental and accessibility gradients". 
# Global Ecology and Biogeography. https://onlinelibrary.wiley.com/doi/full/10.1111/geb.13667
# Codes by Otso Ovaskainen, Mirkka Jones and Gabriela Zuquim
# The codes were modified from November 2020 Hmsc course scripts prepared by Otso Ovaskainen, Jari Oksanen and others.
# Current versions of these training materials are available at https://www.helsinki.fi/en/researchgroups/statistical-ecology/software/hmscCodes
# see readme file for details
##################################################################################################################################################

library(Hmsc)
load("allData.R") #S, X, Y_Adult, Y_Juvenile, Tr, P
#### In the paper (citation above), the juvenile data were not used
rm(Y_Juvenile)

# Check for absent (0) or ubiquitous species (1) in the species data matrices.
range(colMeans(Y_Adult>0))

#Adult Y and Tr data preparation:
# Remove rare species (10 or less observations)
rarespecies = which(colSums(Y_Adult>0)<10)
length(rarespecies)
Y_Adult = Y_Adult[,-rarespecies]

# Species richness distribution (adults only) across sites.
Tr_Adult = droplevels(Tr[row.names(Tr) %in% colnames(Y_Adult),])
Tr_Adult = Tr_Adult[, 1:8]

TrFormula = ~human_food + material + medicine + use_intensity

# prepare environmetal data
X$hab_type2 = as.character(X$hab_type)
X$hab_type2[X$hab_type2=="White sand"] = "Terra firma"
X$hab_type2 = factor(X$hab_type2)
X$log.HAND_50 = log10(X$HAND_50+2)

X = X[, c("Chelsa150_bio6", "Chelsa150_bio15", "Band3_median15", "Band4_median15", "Band7_median15", "logCat", "log.HAND_50", "DistRivlog10", "hab_type2")]

XFormula = ~Chelsa150_bio6 + Chelsa150_bio15 + Band3_median15 + Band4_median15 + Band7_median15 + poly(logCat, 2, raw = TRUE) + log.HAND_50 + DistRivlog10 + hab_type2

# All 430 site codes and coordinates are unique. Coded a spatial random effect.
studyDesign = data.frame(transect = as.factor(S$Transect))

xy = data.frame(x=S$Longitude,y=S$Latitude)
rownames(xy) = studyDesign$transect
rL.transect = HmscRandomLevel(sData = xy, longlat = TRUE)

# Response Y matrix of presence-absence records of the adult palms.
Ypa_AD = 1*(Y_Adult>0)

#####################################################################
m1 = Hmsc(Y=Ypa_AD, XData = X,  XFormula = XFormula,
          TrData = Tr_Adult,TrFormula = TrFormula,
          phyloTree = P,
          distr="probit",
          studyDesign=studyDesign,
          ranLevels={list("transect" = rL.transect)})

models = list(m1)
modelnames = c("presence_absence_ADULT")

save(models,modelnames,file = "unfitted_models.Rdata")

