#S2
library(Hmsc)
load("allData.R") #S, X, Y_Adult, Y_Juvenile, Tr, P

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

#Juvenile Y and Tr data preparation:
rarespecies = which(colSums(Y_Juvenile>0)<5)
length(rarespecies)
Y_Juvenile = Y_Juvenile[,-rarespecies]

# Species richness distribution (juveniles only) across sites.
Tr_Juvenile = droplevels(Tr[row.names(Tr) %in% colnames(Y_Juvenile),])
Tr_Juvenile[, 9:10] = apply(Tr_Juvenile[, 9:10], 2, as.numeric)

TrFormula = ~human_food + material + medicine + use_intensity


#prepare environmetal data
library(vegan)
pcoa = cmdscale(vegdist(1*(Y_Adult>0), binary = TRUE), k = 3, eig = TRUE)
pca_clim = prcomp(X[, 2:19], scale = TRUE)

cortable = round(cor(data.frame(pcoa$points[, 1], pcoa$points[, 2], pcoa$points[, 3], pca_clim$x[, 1:3], X[, -29]), use = "complete.obs"), digits = 3)[, 1:3]# Correlation table sorted by correlations with taxa PCoA axis 1:
head(cortable[order(cortable[, 1]), ]);tail(cortable[order(cortable[, 1]), ])
# Correlation table sorted by correlations with taxa PCoA axis 1:
head(cortable[order(cortable[, 2]), ]);tail(cortable[order(cortable[, 2]), ])
head(cortable[order(cortable[, 3]), ]);tail(cortable[order(cortable[, 3]), ])

pcoa = cmdscale(vegdist(1*(Y_Juvenile>0), binary = TRUE), k = 3, eig = TRUE)
pca_clim = prcomp(X[, 2:19], scale = TRUE)

cortable = round(cor(data.frame(pcoa$points[, 1], pcoa$points[, 2], pcoa$points[, 3], pca_clim$x[, 1:3], X[, -29]), use = "complete.obs"), digits = 3)[, 1:3]# Correlation table sorted by correlations with taxa PCoA axis 1:
head(cortable[order(cortable[, 1]), ]);tail(cortable[order(cortable[, 1]), ])
# Correlation table sorted by correlations with taxa PCoA axis 1:
head(cortable[order(cortable[, 2]), ]);tail(cortable[order(cortable[, 2]), ])
head(cortable[order(cortable[, 3]), ]);tail(cortable[order(cortable[, 3]), ])

# Based on PCoA correlations, selected two climatic variables quantifying variation in temperature and rainfall, three landsat bands, log(cations), HAND, dist to river and habitat type.

X = X[, c("Chelsa150_bio3", "Chelsa150_bio15", "Band2_median15", "Band4_median15", "Band7_median15", "logCat", "HAND_50", "DistRivlog10", "hab_type")]
relevel(X$hab_type, ref = "Terra firma")

XFormula = ~Chelsa150_bio3 + Chelsa150_bio15 + Band2_median15 + Band4_median15 + Band7_median15 + poly(logCat, 2, raw = TRUE) + HAND_50 + DistRivlog10 + hab_type


# All 430 site-seedling codes and coordinates are unique. Coded as a spatial random effect.

studyDesign = data.frame(transect = as.factor(S$Transect))

xy = data.frame(x=S$Longitude,y=S$Latitude)
rownames(xy) = studyDesign$transect
rL.transect = HmscRandomLevel(sData = xy, longlat = TRUE)


Ypa_JUV = 1*(Y_Juvenile>0)
Yabu_JUV = Y_Juvenile
Yabu_JUV[Yabu_JUV==0] = NA
Yabu_JUV=log(Yabu_JUV)

Ypa_AD = 1*(Y_Adult>0)
Yabu_AD = Y_Adult
Yabu_AD[Yabu_AD==0] = NA
Yabu_AD=log(Yabu_AD)

#####################################################################

m1 = Hmsc(Y=Ypa_JUV, XData = X,  XFormula = XFormula,
         TrData = Tr_Juvenile,TrFormula = TrFormula,
         phyloTree = P,
         distr="probit",
         studyDesign=studyDesign,
         ranLevels={list("transect" = rL.transect)})

m2 = Hmsc(Y=Yabu_JUV, YScale = TRUE,
          XData = X,  XFormula = XFormula,
          TrData = Tr_Juvenile,TrFormula = TrFormula,
          phyloTree = P,
          distr="normal",
          studyDesign=studyDesign,
          ranLevels={list("transect" = rL.transect)})

models = list(m1, m2)
modelnames = c("presence_absence_JUV", "abundance_cop_JUV")

save(models,modelnames,file = file.path(ModelDir, "unfitted_modelsJUV"))

#####################################################################

m1 = Hmsc(Y=Ypa_AD, XData = X,  XFormula = XFormula,
          TrData = Tr_Adult,TrFormula = TrFormula,
           phyloTree = P,
           distr="probit",
           studyDesign=studyDesign,
           ranLevels={list("transect" = rL.transect)})

m2 = Hmsc(Y=Yabu_AD, YScale = TRUE,
           XData = X,  XFormula = XFormula,
           TrData = Tr_Adult,TrFormula = TrFormula,
           phyloTree = P,
           distr="normal",
           studyDesign=studyDesign,
           ranLevels={list("transect" = rL.transect)})

models = list(m1, m2)
modelnames = c("presence_absence_ADULT", "abundance_cop_ADULT")

save(models,modelnames,file = "unfitted_models")

