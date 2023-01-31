# HMSC_AmazonPalmUses
Codes associated to the analysis used the paper "The legacy of human use in Amazonian palm communities along environmental and accessibility gradients" by Gabriela Zuquim, Mirkka M. Jones, Otso Ovaskainen, William Trujillo and Henrik Balslev. Global Ecology and Biogeography. In press.

Codes written by Otso Ovaskainen, Mirkka Jones and Gabriela Zuquim
The codes were modified from November 2020 Hmsc course scripts. Current versions of these training materials are available at https://www.helsinki.fi/en/researchgroups/statistical-ecology/software/hmsc

Contact person: gabriela.zuquim@utu.fi 

Raw data was obtained from the sources described below. Please contact the responsible partners for inquiries.

Palm data
We obtained palm species data from inventories in 430 transects (each 500 x 5 m) collected across the study region between 1995 and 2012 (details in Balslev et al., 2019). In each transect, all palm individuals were counted and identified to the species or variety level. The original dataset was quantitative, and classified individuals into seedling, juvenile, adult, and sub-adult categories. Here, we transform the original dataset from abundance to presence-absence data. Moreover, we include only adults and sub-adults, which were combined into a single category. Juveniles and seedlings were excluded from the analysis. Vouchers of collected specimens are deposited at the Aarhus University herbarium (AAU) and in at least one herbarium in the countries of origin of the specimens (COL, QCA, AMAZ, USZ, LPB, INPA). 

Environmental characteristics
Local environmental variables for each transect were obtained from the same published dataset (Balslev et al., 2019). Soil nutrient concentrations (SOIL; the sum of the exchangeable cations Na, Ca, Mg, and K) were measured in the laboratory from topsoil (0-10 cm) samples typically taken at the beginning, middle, and end of each transect. We used a single soil value, the log-transformed exchangeable cation sum averaged over all samples within each transect, to represent its soil fertility in subsequent analyses. Habitat types (HABITAT) were based on field observations, interviews with local guides, and by visual examination of satellite images. The original classification included seven categories condensed into two for analysis. The categories were coded as inundated forests, which include floodplains, backswamps, and restingas (1), and terra-firme forests, which include pre-montane hills, terra-firme and white-sand forests and terraces (2). Anthropogenic soil areas (Amazonian Dark Earth; ADE) were not sampled. 

Hydrological conditions (HAND) were coded as the vertical height of each transect above the nearest drainage (Rennó et al., 2008; Nobre et al., 2011). HAND values were extracted from a publicly available GIS layer at 90 m resolution and based on a 50-pixel drainage area contribution threshold designed to include smaller drainages (http://www.dpi.inpe.br/Ambdata/English/download.php). Nineteen bioclimatic variables (CHELSA) were extracted at 30 arcsec resolution (approximately 1 km) from the Climatologies at High Resolution for the Earth's Land Surface Areas database (Karger et al., 2017). As a measure of site accessibility, we calculated the shortest log-transformed distance from each transect to the nearest river (DISTRIVER) as estimated from the HydroRivers GIS-layer (https://www.hydrosheds.org/page/hydrorivers). Watercourses of orders 1–5 were considered rivers.

Landsat TM/ETM+ reflectance values (LANDSAT) were extracted from a 30 m resolution composite of Landsat images based on over 16,000 relatively cloud‐free acquisitions for 2000–2009 (Van doninck & Tuomisto, 2018). We obtained a single value for each band to represent the whole transect by calculating the median value of all the pixels in a 15 × 15-pixel window centered on the transect. This was done separately for each of the three spectral bands (bands 3, 4 and 7) analyzed in this study.

Regional predictions of palm community attributes across the western Amazon required a regional soil raster as input. To generate a high-quality soil raster, we retrieved field-sampled soil exchangeable base cation (EBC) values from the Harmonized World Soil Database, from databases maintained by the Amazon Research Team of University of Turku (utu.fi/amazon), Brazilian Program from Biodiversity Research (https://ppbio.inpa.gov.br/en), and Balslev et al. (Balslev et al., 2019). The EBC values were also indirectly estimated using EBC-indicator fern species occurrences for the entire region generated following the methods described by Zuquim et al. (Zuquim et al., 2019). The georeferenced soil EBC estimates served as input values in random forest models using median Landsat band 2, 3, 4, 5 and 7 values within a 200 m buffer as predictors. The models were evaluated using random 10-fold cross-validation. The best RF model (R2 = 0.32) was then projected to cover the whole study region using Landsat band values at a 450 m resolution.

Human uses
We recorded each palm species` types and diversity of known human uses, if any, from the World Checklist of useful plants (Diazgranados et al., 2020). The three main human use categories included were food, building material (e.g., housing, thatch), and medicine. For each of these uses, species received a binary classification of used (1) or not used (0). We also assigned a ‘use diversity’ value to each species based on seven types of use, including human food, animal food, material (e.g., housing, thatch), medicine, environmental uses (e.g., ornamental, shading other plants, soil improvers), fuel, genetic, and social uses. These values ranged from 0 to 7. The underlying assumption is that a higher use diversity indicates that a species is, or has historically been more intensively used by humans (Macía et al., 2011). These human use values were used in subsequent analyses to explore whether species with more uses differed in their realized environmental niches or spatial distribution from species with fewer uses.

Phylogeny
We obtained a genus-level phylogenetic tree based on 1000 palm phylogenies (Faurby et al., 2016; Phylogeny_Con_Checklist.nex). The original 1000 trees were reduced to include only the species of interest using the ape (Paradis & Schliep, 2019) and phytools (Revell, 2012) packages in R (R Core Team, 2020). A Maximum Clade Credibility (MCC) tree for the 1000 trees was calculated using TreeAnnotator v2.6.3 (BEAST 2 package (Bouckaert et al., 2014)). All branches with a posterior support below 95% were collapsed using TreeGraph 2 (Stöver & Müller, 2010). All MCC tree nodes showed good support at the genus level but not at the species level. We, therefore, used genus-level data in the analysis.


