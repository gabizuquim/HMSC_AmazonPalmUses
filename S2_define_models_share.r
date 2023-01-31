#####Codes by OO, MJ, GZ####
#####see readme.X for details####

library(Hmsc)
load("allData.R") #S, X, Y_Adult, Y_Juvenile, Tr, P
####For the purpose of the paper CITE, Juvenile data was not used

#### X$HAND_50 has 13 missing values. Exclude these rows from S, X and Y! ####

X$HAND_50 = as.numeric(as.character(X$HAND_50))
missingvalues = which(is.na(X$HAND_50))
S = droplevels(S[-missingvalues,])
X = droplevels(X[-missingvalues,])
Y_Adult = Y_Adult[-missingvalues,]


# Check for absent (0) or ubiquitous species (1) in the species data matrices.
range(colMeans(Y_Adult>0))
min(colSums(Y_Adult>0))
# =0.

#Adult Y and Tr data preparation:
# Remove rare species (5 or less observations)
rarespecies = which(colSums(Y_Adult>0)<5)
length(rarespecies)
Y_Adult = Y_Adult[,-rarespecies]

# Species richness distribution (adults only) across sites.
Tr_Adult = droplevels(Tr[row.names(Tr) %in% colnames(Y_Adult),])
Tr_Adult[, 9:10] = apply(Tr_Adult[, 9:10], 2, as.numeric)

TrFormula = ~human_food + material + medicine + use_intensity

#prepare environmetal data
library(vegan)
pcoa = cmdscale(vegdist(1*(Y_Adult>0), binary = TRUE), k = 3, eig = TRUE)
pca_clim = prcomp(X[, 2:19], scale = TRUE)

# Correlation table sorted by correlations with taxa PCoA axis 1:
cortable = round(cor(data.frame(pcoa$points[, 1], pcoa$points[, 2], pcoa$points[, 3], pca_clim$x[, 1:3], X[, -29]), use = "complete.obs"), digits = 3)[, 1:3]# Correlation table sorted by correlations with taxa PCoA axis 1:
head(cortable[order(cortable[, 1]), ]);tail(cortable[order(cortable[, 1]), ])
head(cortable[order(cortable[, 2]), ]);tail(cortable[order(cortable[, 2]), ])

# Based on PCoA correlations, selected two climatic variables quantifying variation in temperature and rainfall, three landsat bands, log(cations), HAND, dist to river and habitat type.
X = X[, c("Chelsa150_bio3", "Chelsa150_bio15", "Band2_median15", "Band4_median15", "Band7_median15", "logCat", "HAND_50", "DistRivlog10", "hab_type")]
relevel(X$hab_type, ref = "Terra firma")

XFormula = ~Chelsa150_bio3 + Chelsa150_bio15 + Band2_median15 + Band4_median15 + Band7_median15 + poly(logCat, 2, raw = TRUE) + HAND_50 + DistRivlog10 + hab_type

# All 430 site codes and coordinates are unique. Coded as a spatial random effect.
studyDesign = data.frame(transect = as.factor(S$Transect))

xy = data.frame(x=S$Longitude,y=S$Latitude)
rownames(xy) = studyDesign$transect
rL.transect = HmscRandomLevel(sData = xy, longlat = TRUE)

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

