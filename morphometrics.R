load("morphometrics.Rdata")
save.image("morphometrics.Rdata")

# remotes::install_github("phytomosaic/ecole")
# devtools::install_github("chankinonn/GroupStruct")
packages <- c("dplyr", "rstatix", "multcomp", "GroupStruct", "gridExtra",
              "PerformanceAnalytics", "multcomp", "dplyr", "rstatix",
              "RColorBrewer", "ggfortify", "ggplot2", "mclust",
              "clustvarsel", "Hmisc", "ecole", "geomorph")
lapply(packages, library, character.only = TRUE)

# Data formatting: Linear and meristic data ----

pheno_full <- read.csv("data/phenotype/lineal_meristica_all.csv", sep = ",")
View(pheno_full)
pheno_full[2:3] <- lapply(pheno_full[2:3], factor) # coerce multiple columns to factor


# Identify juveniles in each group (optional)
# discard specimens for wich SVL < 70 % of SVL max

# split df by group
# pheno_tax <- split(pheno_full, pheno_full$hyp_tax)

# for(i in names(pheno_tax)){
#   pheno_tax[[i]] <- pheno_tax[[i]] %>% group_by(sex) %>% filter(SVL >= (max(SVL)*0.7))
#   pheno_tax[[i]] <- as.data.frame(pheno_tax[[i]])
# }

# count retained individuals
# x <- c() 
# for(i in names(pheno_tax)){
#   x[i] <- (nrow(pheno_tax[[i]]))
# }
# sum(x)
# setdiff(pheno_full$voucher, bind_rows(pheno_tax)$voucher) # which were the discarded specimens?

# re-join data frames of species
# pheno_full2 <- bind_rows(pheno_tax)

# split in linear and meristic datasets
linear_raw <- pheno_full[, -c(19:26)]
merist_raw <- pheno_full[, -c(4:18, 27)]

# explore variability of some attributes with boxplot
ggplot(linear_raw, aes(x = taxon_sp, y = RH)) + geom_boxplot()
ggplot(linear_raw, aes(x = taxon_sp, y = PDO)) + geom_boxplot()

# Test for correlation -> discard R^2 >= 0.9
# Pearson for linear data and Spearman for meristic data
res_lin <- rcorr(as.matrix(linear_raw[, 4:18]), type = "pearson")
res_mer <- rcorr(as.matrix(merist_raw[, 4:10]), type = "spearman")

# function to format the correlation matrix
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
} # cormat: matrix of the correlation coefficients; pmat: matrix of the correlation p-values

cor_lin <- flattenCorrMatrix(res_lin$r, res_lin$P)
cor_lin[cor_lin$cor >= 0.9, ] # which variables are correlated? -> ADO and PDO -> remove PDO

cor_mer <- flattenCorrMatrix(res_mer$r, res_mer$P)
cor_mer[cor_mer$cor >= 0.9, ] # which variables are correlated? -> none

#### Allometric body size correction ####
# input dataset with no outliers
# for the allom() function to work, remove the "sex" column
# or any other additional columns other than species/population and morphological variables
# we recommend that users apply the multispecies method that can better account for potential differences in average body size among OTUs

linear_allom <- allom(linear_raw[, -c(1, 3, 18:22)], type = "species") # only keep taxonomic designation and uncorrelated variables
# here I removed PDO

# function to normalize the meristic dataset
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
merist_norm <- as.data.frame(lapply(merist_raw[4:10], min_max_norm))

# add columns "sex" and "geom_data", and concatenate with meristic data
# first check the match of both DFs
sum(linear_allom$taxon_sp == linear_raw$taxon_sp) == nrow(linear_allom)
linear_allom$sex <- linear_raw$sex
# relocate column "sex"
linear_allom <- linear_allom %>% relocate(sex, .after = taxon_sp)
# concatenate with meristic data
data_l <- cbind(linear_allom, merist_norm)

#### Test the effects of sex and species with MANOVA ####

summary(manova(lm(as.matrix(data_l[, 3:23]) ~ data_l$sex * data_l$taxon_sp)))
#                                  Df  Pillai approx F num Df den Df Pr(>F)
# data_l$sex                        1 0.17049   6.7826     21    693 <2e-16 ***
# data_l$taxon_sp                  15 1.88383   4.8354    315  10605 <2e-16 ***
# data_l$sex:data_l$taxon_sp       15 0.42918   0.9916    315  10605 0.5312
# Residuals                       713
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# there is a significant effect of sex and species -> create data for each sex
data_sex <- split(data_l, data_l$sex)
data_m <- data_sex$male; data_m[, "sex"] <- NULL
data_f <- data_sex$female; data_f[, "sex"] <- NULL



# Geometric data ----

# read .tps
gm <- readland.tps("data/phenotype/grupoKingii_males.tps",
                   readcurves = TRUE, specID = "imageID")
# use readcurves = TRUE if you digitize curves in TPS

# define semilandmarks
semildm <- c(4, 7:16, 5) # 7 to 16, fixed by landmarks 4 and 5
curve <- matrix(0, nrow = 10, ncol = 3)
for(i in 1:10){
  curve[i, ] <- c(semildm[i], semildm[i + 1], semildm[i + 2])
}

# Procrustes superimposition
proc <- gpagen(gm, curves = curve)

# plot the aligned specimens
plotAllSpecimens(proc$coords)
# check outliers
out <- plotOutliers(proc$coords)
out
omit <- out[1:12] # first 12 are outliers

proc2d <- two.d.array(proc$coords)
proc_red <- proc2d[-omit, ]
# add classifier
classif_gm <- read.csv("data/sheets/classif_gm.csv")
classif_gm <- classif_gm[match(rownames(proc_red), classif_gm$voucher), ]
data_gm <- data.frame(factor(classif_gm$taxon_sp), proc_red)
colnames(data_gm)[1] <- "taxon_sp"


#### PC9 loadings ####
# check which landmarks explain the majority of variability in PC9
# (this landmark showed signals of transgressive evolution)

data <- read.csv("data/phenotype/phenoevol_gm_raw.csv")
pca <- gm.prcomp(data[, -(1:2)])
loadspc9 <- pca$rotation[, 9]
sorted_indices <- order(-abs(loadspc9))
loadspc9[sorted_indices]


#### Deformation grids ####

data_gm_array <- list()
data_gm_array$coords <- readland.tps(
  "data/phenotype/gm_males_red_nooutliers.tps",
  readcurves = TRUE,
  specID = "imageID"
  )
classif <- classif_gm[match(dimnames(data_gm_array$coords)[[3]],
                            classif_gm$voucher), ]
data_gm_array$sp <- classif$taxon_sp 

semildm <- c(4, 7:16, 5) # 7 to 16, fixed by landmarks 4 and 5
curve <- matrix(0, nrow = 10, ncol = 3)
for(i in 1:10){
  curve[i, ] <- c(semildm[i], semildm[i + 1], semildm[i + 2])
}

# split coords based on species
subset_coords <- coords.subset(data_gm_array$coords, data_gm_array$sp)

# here I subset according to the different subgroups defined below for the PCA
data <- list()
data$coords <- abind::abind(subset_coords$baguali,
                            subset_coords$escarchadosi,
                            subset_coords$sarmientoi,
                            subset_coords$tari,
                            subset_coords$spA,
                            along = 3)
proc <- gpagen(data$coords, curves = curve)
PCA <- gm.prcomp(proc$coords)
msh <- mshape(proc$coords)
plotRefToTarget(PCA$shapes$shapes.comp2$min, PCA$shapes$shapes.comp2$max,
                method = "vector", mag = 2,
                gridPars = gridPar(pt.bg = c(rep("red", 6), rep("blue", 10)),
                                   pt.size = 1.5))
par(lwd = 1, col = "black")
box()



# PCA ----

# plot PDF 3 x 4 in. (ldscape)

#### Northern pops: attenboroughi, somuncurae, uptoni ####

# linear data (males)
data_north_m <- filter(data_m,
                       taxon_sp == "somuncurae" |
                       taxon_sp == "uptoni" |
                       taxon_sp == "attenboroughi")
pca_north_m <- prcomp(data_north_m[, -1])
autoplot(pca_north_m, data = data_north_m, colour = "taxon_sp", frame = TRUE) +
  theme_classic(base_family = "LM Sans 10") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")
# linear data (females)
data_north_f <- filter(data_f,
                       taxon_sp == "somuncurae" |
                         taxon_sp == "uptoni" |
                         taxon_sp == "attenboroughi")
pca_north_f <- prcomp(data_north_f[, -1])
autoplot(pca_north_f, data = data_north_f, colour = "taxon_sp", frame = TRUE) +
  theme_classic(base_family = "LM Sans 10") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")
# geometric data (males)
data_north_g <- filter(data_gm,
                       taxon_sp == "somuncurae" |
                         taxon_sp == "uptoni" |
                         taxon_sp == "attenboroughi")
pca_north_g <- data.frame(droplevels(data_north_g$taxon_sp),
                          gm.prcomp(data_north_g[, -1])$x)
colnames(pca_north_g)[1] <- "taxon_sp"
hull_north_g <- pca_north_g[, 1:3] %>%
  group_by(taxon_sp) %>%
  slice(chull(Comp1, Comp2))
pca_north_g %>%
  ggplot(aes(x = Comp1, y = Comp2)) +
  geom_point(aes(x = Comp1, y = Comp2, col = pca_north_g$taxon_sp)) +
  labs(colour = "taxon_sp") +
  geom_polygon(data = hull_north_g, aes(fill = taxon_sp, colour = taxon_sp),
               alpha = 0.3, show.legend = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic(base_family = "LM Sans 10")

#### Southern pops: baguali - sarmientoi - escarchadosi - tari - sp. A ####

# linear data (males)
data_south_m <- filter(data_m,
                       taxon_sp == "sarmientoi" |
                         taxon_sp == "escarchadosi" |
                         taxon_sp == "tari" |
                         taxon_sp == "spA" |
                         taxon_sp == "baguali")
pca_south_m <- prcomp(data_south_m[, -1])
autoplot(pca_south_m, data = data_south_m, colour = "taxon_sp", frame = TRUE) +
  theme_classic(base_family = "LM Sans 10") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")
# linear data (females)
data_south_f <- filter(data_f,
                       taxon_sp == "sarmientoi" |
                         taxon_sp == "escarchadosi" |
                         taxon_sp == "tari" |
                         taxon_sp == "spA" |
                         taxon_sp == "baguali")
pca_south_f <- prcomp(data_south_f[, -1])
autoplot(pca_south_f, data = data_south_f, colour = "taxon_sp", frame = TRUE) +
  theme_classic(base_family = "LM Sans 10") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")
# geometric data (males)
data_south_g <- filter(data_gm,
                       taxon_sp == "sarmientoi" |
                         taxon_sp == "escarchadosi" |
                         taxon_sp == "tari" |
                         taxon_sp == "spA" |
                         taxon_sp == "baguali")
pca_south_g <- data.frame(droplevels(data_south_g$taxon_sp), gm.prcomp(data_south_g[, -1])$x)
colnames(pca_south_g)[1] <- "taxon_sp"
hull_south_g <- pca_south_g[, 1:3] %>%
  group_by(taxon_sp) %>%
  slice(chull(Comp1, Comp2))
pca_south_g %>%
  ggplot(aes(x = Comp1, y = Comp2)) +
  geom_point(aes(x = Comp1, y = Comp2, col = pca_south_g$taxon_sp)) +
  labs(colour = "taxon_sp") +
  geom_polygon(data = hull_south_g, aes(fill = taxon_sp, colour = taxon_sp),
               alpha = 0.3, show.legend = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic(base_family = "LM Sans 10")

#### Western pops: archeforus - scolaroi - zullyae - chacabucoense - tristis ####

# linear data (males)
data_west_m <- filter(data_m,
                      taxon_sp == "archeforus" |
                        taxon_sp == "scolaroi" |
                        taxon_sp == "zullyae" |
                        taxon_sp == "chacabucoense" |
                        taxon_sp == "tristis")
pca_west_m <- prcomp(data_west_m[, -1])
autoplot(pca_west_m, data = data_west_m, colour = "taxon_sp", frame = TRUE) +
  theme_classic(base_family = "LM Sans 10") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")
# linear data (females)
data_west_f <- filter(data_f,
                      taxon_sp == "archeforus" |
                        taxon_sp == "scolaroi" |
                        taxon_sp == "zullyae" |
                        taxon_sp == "chacabucoense" |
                        taxon_sp == "tristis")
pca_west_f <- prcomp(data_west_f[, -1])
autoplot(pca_west_f, data = data_west_f, colour = "taxon_sp", frame = TRUE) +
  theme_classic(base_family = "LM Sans 10") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")
# geometric data (males)
data_west_g <- filter(data_gm,
                      taxon_sp == "archeforus" |
                        taxon_sp == "scolaroi" |
                        taxon_sp == "zullyae" |
                        taxon_sp == "chacabucoense" |
                        taxon_sp == "tristis")
pca_west_g <- data.frame(droplevels(data_west_g$taxon_sp), gm.prcomp(data_west_g[, -1])$x)
colnames(pca_west_g)[1] <- "taxon_sp"
hull_west_g <- pca_west_g[, 1:3] %>%
  group_by(taxon_sp) %>%
  slice(chull(Comp1, Comp2))
pca_west_g %>%
  ggplot(aes(x = Comp1, y = Comp2)) +
  geom_point(aes(x = Comp1, y = Comp2, col = pca_west_g$taxon_sp)) +
  labs(colour = "taxon_sp") +
  geom_polygon(data = hull_west_g, aes(fill = taxon_sp, colour = taxon_sp),
               alpha = 0.3, show.legend = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic(base_family = "LM Sans 10")

#### kingii - "Deseado" ####

# linear data (males)
data_kin_des_m <- filter(data_m, taxon_sp == "kingii" | taxon_sp == "Deseado")
pca_kin_des_m <- prcomp(data_kin_des_m[, -1])
autoplot(pca_kin_des_m, data = data_kin_des_m, colour = "taxon_sp", frame = TRUE) +
  theme_classic(base_family = "lM Sans 10") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")
# linear data (females)
data_kin_des_f <- filter(data_f, taxon_sp == "kingii" | taxon_sp == "Deseado")
pca_kin_des_f <- prcomp(data_kin_des_f[, -1])
autoplot(pca_kin_des_f, data = data_kin_des_f, colour = "taxon_sp", frame = TRUE) +
  theme_classic(base_family = "LM Sans 10") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")
# geometric data (males)
data_kin_des_g <- filter(data_gm, taxon_sp == "kingii" | taxon_sp == "Deseado")
pca_kin_des_g <- data.frame(droplevels(data_kin_des_g$taxon_sp), gm.prcomp(data_kin_des_g[, -1])$x)
colnames(pca_kin_des_g)[1] <- "taxon_sp"
hull_kin_des_g <- pca_kin_des_g[, 1:3] %>%
  group_by(taxon_sp) %>%
  slice(chull(Comp1, Comp2))
pca_kin_des_g %>%
  ggplot(aes(x = Comp1, y = Comp2)) +
  geom_point(aes(x = Comp1, y = Comp2, col = pca_kin_des_g$taxon_sp)) +
  labs(colour = "taxon_sp") +
  geom_polygon(data = hull_kin_des_g, aes(fill = taxon_sp, colour = taxon_sp),
               alpha = 0.3, show.legend = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic(base_family = "LM Sans 10")

#### gallardoi ####
# vs baguali, "Deseado", West cluster (archeforus + scolaroi + zullyae)

# linear data (males)
data_gall_m <- filter(data_m,
                      taxon_sp == "gallardoi" |
                        taxon_sp == "baguali" |
                        taxon_sp == "Deseado"
                      | taxon_sp == "archeforus" |
                        taxon_sp == "scolaroi" |
                        taxon_sp == "zullyae")
levels(data_gall_m$taxon_sp)[c(1, 10, 16)] <- "westclust"
pca_gall_m <- prcomp(data_gall_m[, -1])
autoplot(pca_gall_m, data = data_gall_m, colour = "taxon_sp", frame = TRUE) +
  theme_classic(base_family = "LM Sans 10") + scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")
# linear data (females)
data_gall_f <- filter(data_f,
                      taxon_sp == "gallardoi" |
                        taxon_sp == "baguali" |
                        taxon_sp == "Deseado"
                      | taxon_sp == "archeforus" |
                        taxon_sp == "scolaroi" |
                        taxon_sp == "zullyae")
levels(data_gall_f$taxon_sp)[c(1, 10, 16)] <- "westclust"
pca_gall_f <- prcomp(data_gall_f[, -1])
autoplot(pca_gall_f, data = data_gall_f, colour = "taxon_sp", frame = TRUE) +
  theme_classic(base_family = "LM Sans 10") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")
# geometric data (males)
data_gall_g <- filter(data_gm,
                      taxon_sp == "gallardoi" |
                        taxon_sp == "baguali" |
                        taxon_sp == "Deseado" |
                        taxon_sp == "archeforus" |
                        taxon_sp == "scolaroi" |
                        taxon_sp == "zullyae")
levels(data_gall_g$taxon_sp)[c(1, 10, 16)] <- "westclust"
pca_gall_g <- data.frame(droplevels(data_gall_g$taxon_sp), gm.prcomp(data_gall_g[, -1])$x)
colnames(pca_gall_g)[1] <- "taxon_sp"
hull_gall_g <- pca_gall_g[, 1:3] %>%
  group_by(taxon_sp) %>%
  slice(chull(Comp1, Comp2))
pca_gall_g %>%
  ggplot(aes(x = Comp1, y = Comp2)) +
  geom_point(aes(x = Comp1, y = Comp2, col = pca_gall_g$taxon_sp)) +
  labs(colour = "taxon_sp") +
  geom_polygon(data = hull_gall_g, aes(fill = taxon_sp, colour = taxon_sp),
               alpha = 0.3, show.legend = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic(base_family = "lmsans10")



# PERMANOVA ----
# here I use PC axes explaining up to 90 % of variance
# use summary(pca_north_m) to explore the proportion of variance explained by each PC
# use gm.prcomp(data_north_g[, -1]) to explore the geometric data

#### Northern pops: attenboroughi, somuncurae, uptoni ####

# linear data (males)
pcscores_north_m <- data.frame(taxon_sp = droplevels(data_north_m$taxon_sp), pca_north_m$x)
permanova_pairwise(vegdist(pcscores_north_m[, 2:7], "euclidean"), pcscores_north_m$taxon_sp)
# linear data (females)
pcscores_north_f <- data.frame(taxon_sp = droplevels(data_north_f$taxon_sp), pca_north_f$x)
permanova_pairwise(vegdist(pcscores_north_f[, 2:8], "euclidean"), pcscores_north_f$taxon_sp)
# geometric data
permanova_pairwise(vegdist(pca_north_g[, 2:7], "euclidean"), pca_north_g$taxon_sp)

#### Southern clade: baguali  sarmientoi - escarchadosi - tari - sp. A ####

# linear data (males)
pcscores_south_m <- data.frame(taxon_sp = droplevels(data_south_m$taxon_sp), pca_south_m$x)
permanova_pairwise(vegdist(pcscores_south_m[, 2:8], "euclidean"), pcscores_south_m$taxon_sp)
# linear data (females)
pcscores_south_f <- data.frame(taxon_sp = droplevels(data_south_f$taxon_sp), pca_south_f$x)
permanova_pairwise(vegdist(pcscores_south_f[, 2:7], "euclidean"), pcscores_south_f$taxon_sp)
# geometric data
permanova_pairwise(vegdist(pca_south_g[, 2:9], "euclidean"), pca_south_g$taxon_sp)

### Western clade: archeforus - chacabucoense - scolaroi - tristis - zullyae ####

# linear data (males)
pcscores_west_m <- data.frame(taxon_sp = droplevels(data_west_m$taxon_sp), pca_west_m$x)
permanova_pairwise(vegdist(pcscores_west_m[, 2:8], "euclidean"), pcscores_west_m$taxon_sp)
# linear data (females)
pcscores_west_f <- data.frame(taxon_sp = droplevels(data_west_f$taxon_sp), pca_west_f$x)
permanova_pairwise(vegdist(pcscores_west_f[, 2:8], "euclidean"), pcscores_west_f$taxon_sp)
# geometric data
permanova_pairwise(vegdist(pca_west_g[, 2:10], "euclidean"), pca_west_g$taxon_sp)

#### kingii - "Deseado" clade ####

# linear data (males)
pcscores_kin_des_m <- data.frame(taxon_sp = droplevels(data_kin_des_m$taxon_sp), pca_kin_des_m$x)
permanova_pairwise(vegdist(pcscores_kin_des_m[, 2:8], "euclidean"), pcscores_kin_des_m$taxon_sp)
# linear data (females)
pcscores_kin_des_f <- data.frame(taxon_sp = droplevels(data_kin_des_f$taxon_sp), pca_kin_des_f$x)
permanova_pairwise(vegdist(pcscores_kin_des_f[, 2:8], "euclidean"), pcscores_kin_des_f$taxon_sp)
# geometric data
permanova_pairwise(vegdist(pca_kin_des_g[, 2:10], "euclidean"), pca_kin_des_g$taxon_sp)

#### gallardoi ####
# vs baguali, "Deseado" clade, and West cluster

# linear data (males)
pcscores_gall_m <- data.frame(taxon_sp = droplevels(data_gall_m$taxon_sp), pca_gall_m$x)
permanova_pairwise(vegdist(pcscores_gall_m[, 2:8], "euclidean"), pcscores_gall_m$taxon_sp)
# linear data (females)
pcscores_gall_f <- data.frame(taxon_sp = droplevels(data_gall_f$taxon_sp), pca_gall_f$x)
permanova_pairwise(vegdist(pcscores_gall_f[, 2:8], "euclidean"), pcscores_gall_f$taxon_sp)
# geometric data
permanova_pairwise(vegdist(pca_gall_g[, 2:10], "euclidean"), pca_gall_g$taxon_sp)



# Gaussian mixture models ----
# Fisherian model applicable to (polygenic) continuous traits
# this model assumes gene frequencies in Hardy–Weinberg equilibrium and
# phenotypic variation within species is normally distributed
# We use this approach to evaluate the support of competing taxonomic hypotheses: 

mclust.options(hcUse = "VARS") # check Mclust options and make sure hcUSE = "VARS", otherwise change it
mclust.options()

# check the number of specimens per putative species and for each dataset
sort(table(data_m$taxon_sp))
sort(table(data_f$taxon_sp))
sort(table(data_gm$taxon_sp))

#### Northern pops: attenboroughi, somuncurae, uptoni ####

# Linear data (males)
# variable selection
north_m_varsel <- clustvarsel(pcscores_north_m[, -1],
                              G = 1:3,
                              search = "greedy",
                              direction = "forward",
                              samp = FALSE)
north_m_varsel <- north_m_varsel$subset
# calculate empirical support for models:
# att, som, upt
north_3sp_m <- as_tibble(pcscores_north_m[, -1]) %>%
  dplyr::select(all_of(north_m_varsel)) %>%
  MclustDA(class = pcscores_north_m$taxon_sp, G = 1)
# att, som+upt
north_2sp_hyp_m <- pcscores_north_m$taxon_sp
levels(north_2sp_hyp_m)[c(2, 3)] <- "som_upt"
north_2sp_m <- as_tibble(pcscores_north_m[, -1]) %>%
  dplyr::select(all_of(north_m_varsel)) %>%
  MclustDA(class = north_2sp_hyp_m, G = 1)
# att+som+upt
north_1sp_m <- as_tibble(pcscores_north_m[, -1]) %>%
  dplyr::select(all_of(north_m_varsel)) %>%
  MclustDA(class = rep("north", nrow(pcscores_north_m)), G = 1)
# compare models
sort(c(model3sp = north_3sp_m$bic,  # att, som, upt
       model2sp = north_2sp_m$bic,  # att, som+upt
       model1sp = north_1sp_m$bic)) # att+som+upt
# model3sp model2sp model1sp 
# 149.2844 160.5016 172.1146  

# Linear data (females)
# variable selection
north_f_varsel <- clustvarsel(pcscores_north_f[, -1],
                              G = 1:3,
                              search = "greedy",
                              direction = "forward",
                              samp = FALSE)
north_f_varsel <- north_f_varsel$subset
# calculate empirical support for models:
# att, som, upt
north_3sp_f <- as_tibble(pcscores_north_f[, -1]) %>%
  dplyr::select(all_of(north_f_varsel)) %>%
  MclustDA(class = as.vector(pcscores_north_f$taxon_sp), G = 1)
# att, som+upt
north_2sp_hyp_f <- pcscores_north_f$taxon_sp
levels(north_2sp_hyp_f)[c(2, 3)] <- "som_upt"
north_2sp_f <- as_tibble(pcscores_north_f[, -1]) %>%
  dplyr::select(all_of(north_f_varsel)) %>%
  MclustDA(class = north_2sp_hyp_f, G = 1)
# att+som+upt
north_1sp_f <- as_tibble(pcscores_north_f[, -1]) %>%
  dplyr::select(all_of(north_f_varsel)) %>%
  MclustDA(class = rep("north", nrow(pcscores_north_f)), G = 1)
# compare models
sort(c(model3sp = north_3sp_f$bic,  # att, som, upt
       model2sp = north_2sp_f$bic,  # att, som+upt
       model1sp = north_1sp_f$bic)) # att+som+upt
# model3sp model2sp model1sp 
# 221.7312 232.3423 243.1785 

# Geometric data (males)
# variable selection
north_g_varsel <- clustvarsel(pca_north_g[, -1],
                              G = 1:3,
                              search = "greedy",
                              direction = "forward",
                              samp = FALSE)
north_g_varsel <- north_g_varsel$subset
# calculate empirical support for models:
# att, som, upt
north_3sp_g <- as_tibble(pca_north_g[, -1]) %>%
  dplyr::select(all_of(north_g_varsel)) %>%
  MclustDA(class = pca_north_g$taxon_sp, G = 1)
# att, som+upt
north_2sp_hyp_g <- pca_north_g$taxon_sp
levels(north_2sp_hyp_g)[c(2, 3)] <- "som_upt"
north_2sp_g <- as_tibble(pca_north_g[, -1]) %>%
  dplyr::select(all_of(north_g_varsel)) %>%
  MclustDA(class = north_2sp_hyp_g, G = 1)
# att+som+upt
north_1sp_g <- as_tibble(pca_north_g[, -1]) %>%
  dplyr::select(all_of(north_g_varsel)) %>%
  MclustDA(class = rep("north", nrow(pca_north_g)), G = 1)
# compare models
sort(c(model3sp = north_3sp_g$bic,  # att, som, upt
       model2sp = north_2sp_g$bic,  # att, som+upt
       model1sp = north_1sp_g$bic)) # att+som+upt
# model3sp model2sp model1sp 
# 346.7372 363.9858 373.0306

#### Southern pops: baguali - escarchadosi - sarmientoi - tari - sp. A ####

# Linear data (males)
# variable selection
south_m_varsel <- clustvarsel(pcscores_south_m[, -1],
                              G = 1:5,
                              search = "greedy",
                              direction = "forward",
                              samp = FALSE)
south_m_varsel <- south_m_varsel$subset
# calculate empirical support for models:
# bag, esc, sar, tar, spA
south_5sp_m <- pcscores_south_m[, -1] %>%
  dplyr::select(all_of(south_m_varsel)) %>%
  MclustDA(class = pcscores_south_m$taxon_sp, G = 1)
# bag, esc+spA, sar, tar
south_4sp_hyp_m <- pcscores_south_m$taxon_sp
levels(south_4sp_hyp_m)[c(2, 4)] <- "esc_spA"
south_4sp_m <- pcscores_south_m[, -1] %>%
  dplyr::select(all_of(south_m_varsel)) %>%
  MclustDA(class = as.vector(south_4sp_hyp_m), G = 1)
# bag, esc+spA+sar, tar
south_3Asp_hyp_m <- south_4sp_hyp_m
levels(south_3Asp_hyp_m)[c(2, 3)] <- "esc_spA_sarm"
south_3Asp_m <- pcscores_south_m[, -1] %>%
  dplyr::select(all_of(south_m_varsel)) %>%
  MclustDA(class = south_3Asp_hyp_m, G = 1)
# bag, esc+spA+tar, sar
south_3Bsp_hyp_m <- south_4sp_hyp_m
levels(south_3Bsp_hyp_m)[c(2, 4)] <- "esc_spA_tar"
south_3Bsp_m <- pcscores_south_m[, -1] %>%
  dplyr::select(all_of(south_m_varsel)) %>%
  MclustDA(class = south_3Bsp_hyp_m, G = 1)
# bag, esc+spA+tar+sar
south_2sp_hyp_m <- south_3Bsp_hyp_m
levels(south_2sp_hyp_m)[c(2, 3)] <- "esc_spA_tar_sar"
south_2sp_m <- pcscores_south_m[, -1] %>%
  dplyr::select(all_of(south_m_varsel)) %>%
  MclustDA(class = south_2sp_hyp_m, G = 1)
# bag+esc+spA+tar+sar
south_1sp_m <- pcscores_south_m[, -1] %>%
  dplyr::select(all_of(south_m_varsel)) %>%
  MclustDA(class = rep("south", nrow(pcscores_south_m)), G = 1)
# compare models
sort(c(model5sp = south_5sp_m$bic,   # bag, esc, spA, sar, tar
       model4sp = south_4sp_m$bic,   # bag, esc+spA, sar, tar
       model3Asp = south_3Asp_m$bic, # bag, esc+spA+sar, tar
       model3Bsp = south_3Bsp_m$bic, # bag, esc+spA+tar, sar
       model2sp = south_2sp_m$bic,   # bag, esc+spA+sar+tar
       model1sp = south_1sp_m$bic))  # bag+esc+spA+sar+tar
# model5sp  model4sp model3Asp model3Bsp  model1sp  model2sp 
# 83.25594 107.35288 121.81861 144.59041 157.98142 158.77267  

# Linear data (females)
# variable selection
pcscores_south_f <- data.frame(taxon_sp = droplevels(data_south_f$taxon_sp), pca_south_f$x)
south_f_varsel <- clustvarsel(pcscores_south_f[, -1],
                              G = 1:5,
                              search = "greedy",
                              direction = "forward",
                              samp = FALSE)
south_f_varsel <- south_f_varsel$subset
# calculate empirical support for models:
# bag, esc, sar, tar, spA
south_5sp_f <- pcscores_south_f[, -1] %>%
  dplyr::select(all_of(south_f_varsel)) %>%
  MclustDA(class = pcscores_south_f$taxon_sp, G = 1)
# bag, esc+spA, sar, tar
south_4sp_hyp_f <- pcscores_south_f$taxon_sp
levels(south_4sp_hyp_f)[c(2, 4)] <- "esc_spA"
south_4sp_f <- pcscores_south_f[, -1] %>%
  dplyr::select(all_of(south_f_varsel)) %>%
  MclustDA(class = south_4sp_hyp_f, G = 1)
# bag, esc+spA+sar, tar
south_3Asp_hyp_f <- south_4sp_hyp_f
levels(south_3Asp_hyp_f)[c(2, 3)] <- "esc_spA_sarm"
south_3Asp_f <- pcscores_south_f[, -1] %>%
  dplyr::select(all_of(south_f_varsel)) %>%
  MclustDA(class = south_3Asp_hyp_f, G = 1)
# bag, esc+spA+tar, sar
south_3Bsp_hyp_f <- south_4sp_hyp_f
levels(south_3Bsp_hyp_f)[c(2, 4)] <- "esc_spA_tar"
south_3Bsp_f <- pcscores_south_f[, -1] %>%
  dplyr::select(all_of(south_f_varsel)) %>%
  MclustDA(class = south_3Bsp_hyp_f, G = 1)
# bag, esc+spA+tar+sar
south_2sp_hyp_f <- south_3Bsp_hyp_f
levels(south_2sp_hyp_f)[c(2, 3)] <- "esc_spA_tar_sar"
south_2sp_f <- pcscores_south_f[, -1] %>%
  dplyr::select(all_of(south_f_varsel)) %>%
  MclustDA(class = south_2sp_hyp_f, G = 1)
# bag+esc+spA+tar+sar
south_1sp_f <- pcscores_south_f[, -1] %>%
  dplyr::select(all_of(south_f_varsel)) %>%
  MclustDA(class = rep("south", nrow(pcscores_south_f)), G = 1)
# compare models
sort(c(model5sp = south_5sp_f$bic,   # bag, esc, spA, sar, tar
       model4sp = south_4sp_f$bic,   # bag, esc+spA, sar, tar
       model3Asp = south_3Asp_f$bic, # bag, esc+spA+sar, tar
       model3Bsp = south_3Bsp_f$bic, # bag, esc+spA+tar, sar
       model2sp = south_2sp_f$bic,   # bag, esc+spA+sar+tar
       model1sp = south_1sp_f$bic))  # bag+esc+spA+sar+tar
# model5sp  model4sp model3Bsp model3Asp  model2sp  model1sp 
# 492.2671  520.5079  567.8779  596.0320  643.8700  681.5301

# Geometric data (males)
# variable selection
south_g_varsel <- clustvarsel(pca_south_g[, -1],
                              G = 1:5,
                              search = "greedy",
                              direction = "forward",
                              samp = FALSE)
south_g_varsel <- south_g_varsel$subset
# calculate empirical support for models:
# bag, esc, sar, tar, spA
south_5sp_g <- pca_south_g[, -1] %>%
  dplyr::select(all_of(south_g_varsel)) %>%
  MclustDA(class = pca_south_g$taxon_sp, G = 1)
# bag, esc+spA, sar, tar
south_4sp_hyp_g <- pca_south_g$taxon_sp
levels(south_4sp_hyp_g)[c(2, 4)] <- "esc_spA"
south_4sp_g <- pca_south_g[, -1] %>%
  dplyr::select(all_of(south_g_varsel)) %>%
  MclustDA(class = south_4sp_hyp_g, G = 1)
# bag, esc+spA+sar, tar
south_3Asp_hyp_g <- south_4sp_hyp_g
levels(south_3Asp_hyp_g)[c(2, 3)] <- "esc_spA_sarm"
south_3Asp_g <- pca_south_g[, -1] %>%
  dplyr::select(all_of(south_g_varsel)) %>%
  MclustDA(class = south_3Asp_hyp_g, G = 1)
# bag, esc+spA+tar, sar
south_3Bsp_hyp_g <- south_4sp_hyp_g
levels(south_3Bsp_hyp_g)[c(2, 4)] <- "esc_spA_tar"
south_3Bsp_g <- pca_south_g[, -1] %>%
  dplyr::select(all_of(south_g_varsel)) %>%
  MclustDA(class = south_3Bsp_hyp_g, G = 1)
# bag, esc+spA+tar+sar
south_2sp_hyp_g <- south_3Bsp_hyp_g
levels(south_2sp_hyp_g)[c(2, 3)] <- "esc_spA_tar_sar"
south_2sp_g <- pca_south_g[, -1] %>%
  dplyr::select(all_of(south_g_varsel)) %>%
  MclustDA(class = south_2sp_hyp_g, G = 1)
# bag+esc+spA+tar+sar
south_1sp_g <- pca_south_g[, -1] %>%
  dplyr::select(all_of(south_g_varsel)) %>%
  MclustDA(class = rep("south", nrow(pca_south_g)), G = 1)
# compare models
sort(c(model5sp = south_5sp_g$bic,   # bag, esc, spA, sar, tar
       model4sp = south_4sp_g$bic,   # bag, esc+spA, sar, tar
       model3Asp = south_3Asp_g$bic, # bag, esc+spA+sar, tar
       model3Bsp = south_3Bsp_g$bic, # bag, esc+spA+tar, sar
       model2sp = south_2sp_g$bic,   # bag, esc+spA+sar+tar
       model1sp = south_1sp_g$bic))  # bag+esc+spA+sar+tar
# model5sp  model4sp model3Bsp model3Asp  model2sp  model1sp 
# 1302.963  1321.993  1338.187  1338.782  1354.245  1374.779

#### Western pops: archeforus - chacabucoense - scolaroi - tristis - zullyae ####

# Linear data (males)
# variable selection
west_m_varsel <- clustvarsel(pcscores_west_m[, -1],
                             G = 1:5,
                             search = "greedy",
                             direction = "forward",
                             emModels2 = mclust.options("emModelNames"),
                             samp = FALSE)
west_m_varsel <- west_m_varsel$subset
# calculate empirical support for models:
# arc, cha, sco, tri, zul
west_5sp_m <- as_tibble(pcscores_west_m[, -1]) %>%
  dplyr::select(all_of(west_m_varsel)) %>%
  MclustDA(class = pcscores_west_m$taxon_sp, G = 1)
# arc, cha, tri, sco+zul
west_4sp_hyp_m <- pcscores_west_m$taxon_sp
levels(west_4sp_hyp_m)[c(3, 5)] <- "sco_zul"
west_4sp_m <- as_tibble(pcscores_west_m[, -1]) %>%
  dplyr::select(all_of(west_m_varsel)) %>%
  MclustDA(class = west_4sp_hyp_m, G = 1)
# arc+sco+zul, cha, tri
west_3Asp_hyp_m <- west_4sp_hyp_m
levels(west_3Asp_hyp_m)[c(1, 3)] <- "arc_sco_zul"
west_3Asp_m <- as_tibble(pcscores_west_m[, -1]) %>%
  dplyr::select(all_of(west_m_varsel)) %>%
  MclustDA(class = west_3Asp_hyp_m, G = 1)
# arc+tri, cha, sco+zul
west_3Bsp_hyp_m <- west_4sp_hyp_m
levels(west_3Bsp_hyp_m)[c(1, 4)] <- "arc_tri"
west_3Bsp_m <- as_tibble(pcscores_west_m[, -1]) %>%
  dplyr::select(all_of(west_m_varsel)) %>%
  MclustDA(class = west_3Bsp_hyp_m, G = 1)
# arc+sco+tri+zul, cha
west_2Asp_hyp_m <- west_3Bsp_hyp_m
levels(west_2Asp_hyp_m)[c(1, 3)] <- "arc_tri_sco_zul"
west_2Asp_m <- as_tibble(pcscores_west_m[, -1]) %>%
  dplyr::select(all_of(west_m_varsel)) %>%
  MclustDA(class = west_2Asp_hyp_m, G = 1)
# arc+cha+sco+zul, tri
west_2Bsp_hyp_m <- west_3Asp_hyp_m
levels(west_2Bsp_hyp_m)[c(1, 2)] <- "arc_cha_sco_zul"
west_2Bsp_m <- as_tibble(pcscores_west_m[, -1]) %>%
  dplyr::select(all_of(west_m_varsel)) %>%
  MclustDA(class = as.vector(west_2Bsp_hyp_m), G = 1)
# arc+tri+sco+zul
west_1sp_m <- as_tibble(pcscores_west_m[, -1]) %>%
  dplyr::select(all_of(west_m_varsel)) %>%
  MclustDA(class = rep("west", nrow(pcscores_west_m)), G = 1)
# compare models
sort(c(model5sp = west_5sp_m$bic,   # arc, cha, sco, tri, zul
       model4sp = west_4sp_m$bic,   # arc, cha, sco+zul, tri
       model3Asp = west_3Asp_m$bic, # arc+sco+zul, cha, tri
       model3Bsp = west_3Bsp_m$bic, # arc+tri, cha, sco+zul
       model2Asp = west_2Asp_m$bic, # arc+tri+sco+zul, cha
       model2Bsp = west_2Bsp_m$bic, # arc+cha+sco+zul, tri
       model1sp = west_1sp_m$bic))  # arc+cha+tri+sco+zul
# model5sp  model4sp model3Bsp model3Asp model2Bsp model2Asp  model1sp 
# 850.9655  921.6919  983.9851  988.4888 1048.2480 1049.0858 1113.8656

# Linear data (females)
# variable selection
west_f_varsel <- clustvarsel(pcscores_west_f[, -1],
                             G = 1:5,
                             search = "greedy",
                             direction = "forward",
                             emModels2 = mclust.options("emModelNames"),
                             samp = FALSE)
west_f_varsel <- west_f_varsel$subset
# calculate empirical support for models:
# arc, cha, sco, tri, zul
west_5sp_f <- as_tibble(pcscores_west_f[, -1]) %>%
  dplyr::select(all_of(west_f_varsel)) %>%
  MclustDA(class = pcscores_west_f$taxon_sp, G = 1)
# arc, cha, tri, sco+zul
west_4sp_hyp_f <- pcscores_west_f$taxon_sp
levels(west_4sp_hyp_f)[c(3, 5)] <- "sco_zul"
west_4sp_f <- as_tibble(pcscores_west_f[, -1]) %>%
  dplyr::select(all_of(west_f_varsel)) %>%
  MclustDA(class = west_4sp_hyp_f, G = 1)
# arc+sco+zul, cha, tri
west_3Asp_hyp_f <- west_4sp_hyp_f
levels(west_3Asp_hyp_f)[c(1, 3)] <- "arc_sco_zul"
west_3Asp_f <- as_tibble(pcscores_west_f[, -1]) %>%
  dplyr::select(all_of(west_f_varsel)) %>%
  MclustDA(class = west_3Asp_hyp_f, G = 1)
# arc+tri, cha, sco+zul
west_3Bsp_hyp_f <- west_4sp_hyp_f
levels(west_3Bsp_hyp_f)[c(1, 4)] <- "arc_tri"
west_3Bsp_f <- as_tibble(pcscores_west_f[, -1]) %>%
  dplyr::select(all_of(west_f_varsel)) %>%
  MclustDA(class = west_3Bsp_hyp_f, G = 1)
# arc+sco+tri+zul, cha
west_2Asp_hyp_f <- west_3Bsp_hyp_f
levels(west_2Asp_hyp_f)[c(1, 3)] <- "arc_tri_sco_zul"
west_2Asp_f <- as_tibble(pcscores_west_f[, -1]) %>%
  dplyr::select(all_of(west_f_varsel)) %>%
  MclustDA(class = west_2Asp_hyp_f, G = 1)
# arc+cha+sco+zul, tri
west_2Bsp_hyp_f <- west_3Asp_hyp_f
levels(west_2Bsp_hyp_f)[c(1, 2)] <- "arc_cha_sco_zul"
west_2Bsp_f <- as_tibble(pcscores_west_f[, -1]) %>%
  dplyr::select(all_of(west_f_varsel)) %>%
  MclustDA(class = west_2Bsp_hyp_f, G = 1)
# arc+tri+sco+zul
west_1sp_f <- as_tibble(pcscores_west_f[, -1]) %>%
  dplyr::select(all_of(west_f_varsel)) %>%
  MclustDA(class = rep("west", nrow(pcscores_west_f)), G = 1)
# compare models
sort(c(model5sp = west_5sp_f$bic,   # arc, cha, sco, tri, zul
       model4sp = west_4sp_f$bic,   # arc, cha, sco+zul, tri
       model3Asp = west_3Asp_f$bic, # arc+sco+zul, cha, tri
       model3Bsp = west_3Bsp_f$bic, # arc+tri, cha, sco+zul
       model2Asp = west_2Asp_f$bic, # arc+tri+sco+zul, cha
       model2Bsp = west_2Bsp_f$bic, # arc+cha+sco+zul, tri
       model1sp = west_1sp_f$bic))  # arc+cha+tri+sco+zul
# model5sp  model4sp model3Bsp model3Asp model2Bsp model2Asp  model1sp 
# 4.034546 18.573399 29.757271 30.744446 39.463391 42.382077 53.899699

# Geometric data (males)
# variable selection
west_g_varsel <- clustvarsel(pca_west_g[, -1],
                             G = 1:5,
                             search = "greedy",
                             direction = "forward",
                             emModels2 = mclust.options("emModelNames"),
                             samp = FALSE)
west_g_varsel <- west_g_varsel$subset
# calculate empirical support for models:
# arc, cha, sco, tri, zul
west_5sp_g <- as_tibble(pca_west_g[, -1]) %>%
  dplyr::select(all_of(west_g_varsel)) %>%
  MclustDA(class = pca_west_g$taxon_sp, G = 1)
# arc, cha, tri, sco+zul
west_4sp_hyp_g <- pca_west_g$taxon_sp
levels(west_4sp_hyp_g)[c(3, 5)] <- "sco_zul"
west_4sp_g <- as_tibble(pca_west_g[, -1]) %>%
  dplyr::select(all_of(west_g_varsel)) %>%
  MclustDA(class = west_4sp_hyp_g, G = 1)
# arc+sco+zul, cha, tri
west_3Asp_hyp_g <- west_4sp_hyp_g
levels(west_3Asp_hyp_g)[c(1, 3)] <- "arc_sco_zul"
west_3Asp_g <- as_tibble(pca_west_g[, -1]) %>%
  dplyr::select(all_of(west_g_varsel)) %>%
  MclustDA(class = west_3Asp_hyp_g, G = 1)
# arc+tri, cha, sco+zul
west_3Bsp_hyp_g <- west_4sp_hyp_g
levels(west_3Bsp_hyp_g)[c(1, 4)] <- "arc_tri"
west_3Bsp_g <- as_tibble(pca_west_g[, -1]) %>%
  dplyr::select(all_of(west_g_varsel)) %>%
  MclustDA(class = west_3Bsp_hyp_g, G = 1)
# arc+sco+tri+zul, cha
west_2Asp_hyp_g <- west_3Bsp_hyp_g
levels(west_2Asp_hyp_g)[c(1, 3)] <- "arc_tri_sco_zul"
west_2Asp_g <- as_tibble(pca_west_g[, -1]) %>%
  dplyr::select(all_of(west_g_varsel)) %>%
  MclustDA(class = west_2Asp_hyp_g, G = 1)
# arc+cha+sco+zul, tri
west_2Bsp_hyp_g <- west_3Asp_hyp_g
levels(west_2Bsp_hyp_g)[c(1, 3)] <- "arc_cha_sco_zul"
west_2Bsp_g <- as_tibble(pca_west_g[, -1]) %>%
  dplyr::select(all_of(west_g_varsel)) %>%
  MclustDA(class = west_2Bsp_hyp_g, G = 1)
# arc+tri+sco+zul
west_1sp_g <- as_tibble(pca_west_g[, -1]) %>%
  dplyr::select(all_of(west_g_varsel)) %>%
  MclustDA(class = rep("west", nrow(pca_west_g)), G = 1)
# compare models
sort(c(model5sp = west_5sp_g$bic,   # arc, cha, sco, tri, zul
       model4sp = west_4sp_g$bic,   # arc, cha, sco+zul, tri
       model3Asp = west_3Asp_g$bic, # arc+sco+zul, cha, tri
       model3Bsp = west_3Bsp_g$bic, # arc+tri, cha, sco+zul
       model2Asp = west_2Asp_g$bic, # arc+tri+sco+zul, cha
       model2Bsp = west_2Bsp_g$bic, # arc+cha+sco+zul, tri
       model1sp = west_1sp_g$bic))  # arc+cha+tri+sco+zul
# model5sp  model4sp model3Asp model3Bsp model2Asp model2Bsp  model1sp 
# 2211.887  2223.233  2230.849  2240.145  2241.105  2249.362  2260.170

#### kingii vs "Deseado" ####

# Linear data (males)
data_kin_des_m <- filter(data_m, taxon_sp == "kingii" | taxon_sp == "Deseado")
pcscores_kin_des_m <- data.frame(taxon_sp = droplevels(data_kin_des_m$taxon_sp),
                                 prcomp(data_kin_des_m[, -1])$x)
# variable selection
kin_des_m_varsel <- clustvarsel(pcscores_kin_des_m[, -1],
                                G = 1:2, search = "greedy",
                                direction = "forward",
                                samp = FALSE)
kin_des_m_varsel <- kin_des_m_varsel$subset
# calculate empirical support for models:
# kin, Des
kin_des_2sp_m <- as_tibble(pcscores_kin_des_m[, -1]) %>%
  dplyr::select(all_of(kin_des_m_varsel)) %>%
  MclustDA(class = pcscores_kin_des_m$taxon_sp, G = 1)
# kin+Des
kin_des_1sp_m <- as_tibble(pcscores_kin_des_m[, -1]) %>%
  dplyr::select(all_of(kin_des_m_varsel)) %>%
  MclustDA(class = rep("kin_des", nrow(pcscores_kin_des_m)), G = 1)
# compare models
sort(c(model2sp = kin_des_2sp_m$bic,  # kin, Des
       model1sp = kin_des_1sp_m$bic)) # kin+Des
# model2sp model1sp 
# 7159.524 7312.339

# Linear data (females)
data_kin_des_f <- filter(data_f, taxon_sp == "kingii" | taxon_sp == "Deseado")
pcscores_kin_des_f <- data.frame(taxon_sp = droplevels(data_kin_des_f$taxon_sp),
                                 prcomp(data_kin_des_f[, -1])$x)
# variable selection
kin_des_f_varsel <- clustvarsel(pcscores_kin_des_f[, -1],
                                G = 1:2,
                                search = "greedy",
                                direction = "forward",
                                samp = FALSE)
kin_des_f_varsel <- kin_des_f_varsel$subset
# calculate empirical support for models:
# kin, Des
kin_des_2sp_f <- as_tibble(pcscores_kin_des_f[, -1]) %>%
  dplyr::select(all_of(kin_des_f_varsel)) %>%
  MclustDA(class = pcscores_kin_des_f$taxon_sp, G = 1)
# kin+Des
kin_des_1sp_f <- as_tibble(pcscores_kin_des_f[, -1]) %>%
  dplyr::select(all_of(kin_des_f_varsel)) %>%
  MclustDA(class = rep("kin_des", nrow(pcscores_kin_des_f)), G = 1)
# compare models
sort(c(model2sp = kin_des_2sp_f$bic,  # kin, Des
       model1sp = kin_des_1sp_f$bic)) # kin+Des
# model2sp model1sp 
# 3550.192 3656.685

# Geometric data (males)
data_kin_des_g <- filter(data_gm, taxon_sp == "kingii" | taxon_sp == "Deseado")
pca_kin_des_g <- data.frame(droplevels(data_kin_des_g$taxon_sp),
                            gm.prcomp(data_kin_des_g[, -1])$x)
colnames(pca_kin_des_g)[1] <- "taxon_sp"
# variable selection
kin_des_g_varsel <- clustvarsel(pca_kin_des_g[, -1],
                                G = 1:2,
                                search = "greedy",
                                direction = "forward",
                                samp = FALSE)
kin_des_g_varsel <- kin_des_g_varsel$subset
# calculate empirical support for models:
# kin, Des
kin_des_2sp_g <- as_tibble(pca_kin_des_g[, -1]) %>%
  dplyr::select(all_of(kin_des_g_varsel)) %>%
  MclustDA(class = pca_kin_des_g$taxon_sp, G = 1)
# kin+Des
kin_des_1sp_g <- as_tibble(pca_kin_des_g[, -1]) %>%
  dplyr::select(all_of(kin_des_g_varsel)) %>%
  MclustDA(class = rep("kin_des", nrow(pca_kin_des_g)), G = 1)
# compare models
sort(c(model2sp = kin_des_2sp_g$bic,  # kin, Des
       model1sp = kin_des_1sp_g$bic)) # kin+Des
# model2sp model1sp 
# 19427.67 19496.28

#### gallardoi ####
#### vs baguali

# Linear data (males)
data_gal_bag_m <- filter(data_m, taxon_sp == "gallardoi" | taxon_sp == "baguali")
pcscores_gal_bag_m <- data.frame(taxon_sp = droplevels(data_gal_bag_m$taxon_sp),
                                 prcomp(data_gal_bag_m[, -1])$x)
# variable selection
gal_bag_m_varsel <- clustvarsel(pcscores_gal_bag_m[, -1],
                                G = 1:2,
                                search = "greedy",
                                direction = "forward", samp = FALSE)
gal_bag_m_varsel <- gal_bag_m_varsel$subset
# calculate empirical support for models:
# gal, bag
gal_bag_2sp_m <- as_tibble(pcscores_gal_bag_m[, -1]) %>%
  dplyr::select(all_of(gal_bag_m_varsel)) %>%
  MclustDA(class = pcscores_gal_bag_m$taxon_sp, G = 1)
# gal+bag
gal_bag_1sp_m <- as_tibble(pcscores_gal_bag_m[, -1]) %>%
  dplyr::select(all_of(gal_bag_m_varsel)) %>%
  MclustDA(class = rep("gal_bag", nrow(pcscores_gal_bag_m)), G = 1)
# compare models
sort(c(model2sp = gal_bag_2sp_m$bic,  # gal, bag
       model1sp = gal_bag_1sp_m$bic)) # gal+bag
# model2sp  model1sp 
# 97.85888 122.40019

# Linear data (females)
data_gal_bag_f <- filter(data_f, taxon_sp == "gallardoi" | taxon_sp == "baguali")
pcscores_gal_bag_f <- data.frame(taxon_sp = droplevels(data_gal_bag_f$taxon_sp),
                                 prcomp(data_gal_bag_f[, -1])$x)
# variable selection
gal_bag_f_varsel <- clustvarsel(pcscores_gal_bag_f[, -1],
                                G = 1:2,
                                search = "greedy",
                                direction = "forward",
                                samp = FALSE)
gal_bag_f_varsel <- gal_bag_f_varsel$subset
# calculate empirical support for models:
# gal, bag
gal_bag_2sp_f <- as_tibble(pcscores_gal_bag_f[, -1]) %>%
  dplyr::select(all_of(gal_bag_f_varsel)) %>%
  MclustDA(class = pcscores_gal_bag_f$taxon_sp, G = 1)
# gal+bag
gal_bag_1sp_f <- as_tibble(pcscores_gal_bag_f[, -1]) %>%
  dplyr::select(all_of(gal_bag_f_varsel)) %>%
  MclustDA(class = rep("gal_bag", nrow(pcscores_gal_bag_f)), G = 1)
# compare models
sort(c(model2sp = gal_bag_2sp_f$bic,  # gal, bag
       model1sp = gal_bag_1sp_f$bic)) # gal+bag
# model2sp model1sp 
# 179.5729 192.0983

# Geometric data (males)
data_gal_bag_g <- filter(data_gm, taxon_sp == "gallardoi" | taxon_sp == "baguali")
pca_gal_bag_g <- data.frame(droplevels(data_gal_bag_g$taxon_sp),
                            gm.prcomp(data_gal_bag_g[, -1])$x)
colnames(pca_gal_bag_g)[1] <- "taxon_sp"
# variable selection
gal_bag_g_varsel <- clustvarsel(pca_gal_bag_g[, -1],
                                G = 1:2,
                                search = "greedy",
                                direction = "forward",
                                samp = FALSE)
gal_bag_g_varsel <- gal_bag_g_varsel$subset  
# calculate empirical support for models:
# gal, bag
gal_bag_2sp_g <- as_tibble(pca_gal_bag_g[, -1]) %>%
  dplyr::select(all_of(gal_bag_g_varsel)) %>%
  MclustDA(class = pca_gal_bag_g$taxon_sp, G = 1)
# gal+bag
gal_bag_1sp_g <- as_tibble(pca_gal_bag_g[, -1]) %>%
  dplyr::select(all_of(gal_bag_g_varsel)) %>%
  MclustDA(class = rep("gal_bag", nrow(pca_gal_bag_g)), G = 1)
# compare models
sort(c(model2sp = gal_bag_2sp_g$bic,  # gal, bag
       model1sp = gal_bag_1sp_g$bic)) # gal+bag
# model2sp model1sp 
# 638.9410 653.6662

#### vs "Deseado" clade

# Linear data (males)
data_gal_des_m <- filter(data_m, taxon_sp == "gallardoi" | taxon_sp == "Deseado")
pcscores_gal_des_m <- data.frame(taxon_sp = droplevels(data_gal_des_m$taxon_sp),
                                 prcomp(data_gal_des_m[, -1])$x)
# variable selection
gal_des_m_varsel <- clustvarsel(pcscores_gal_des_m[, -1],
                                G = 1:2,
                                search = "greedy",
                                direction = "forward",
                                samp = FALSE)
gal_des_m_varsel <- gal_des_m_varsel$subset
# calculate empirical support for models:
# gal, Des
gal_des_2sp_m <- as_tibble(pcscores_gal_des_m[, -1]) %>%
  dplyr::select(all_of(gal_des_m_varsel)) %>%
  MclustDA(class = pcscores_gal_des_m$taxon_sp, G = 1)
# gal+Des
gal_des_1sp_m <- as_tibble(pcscores_gal_des_m[, -1]) %>%
  dplyr::select(all_of(gal_des_m_varsel)) %>%
  MclustDA(class = rep("gal_des", nrow(pcscores_gal_des_m)), G = 1)
# compare models
sort(c(model2sp = gal_des_2sp_m$bic,  # gal, Des
       model1sp = gal_des_1sp_m$bic)) # gal+Des
# model2sp model1sp 
# 908.2652 953.0603

# Linear data (females)
data_gal_des_f <- filter(data_f, taxon_sp == "gallardoi" | taxon_sp == "Deseado")
pcscores_gal_des_f <- data.frame(taxon_sp = droplevels(data_gal_des_f$taxon_sp),
                                 prcomp(data_gal_des_f[, -1])$x)
# variable selection
gal_des_f_varsel <- clustvarsel(pcscores_gal_des_f[, -1],
                                G = 1:2,
                                search = "greedy",
                                direction = "forward",
                                samp = FALSE)
gal_des_f_varsel <- gal_des_f_varsel$subset
# calculate empirical support for models:
# gal, Des
gal_des_2sp_f <- as_tibble(pcscores_gal_des_f[, -1]) %>%
  dplyr::select(all_of(gal_des_f_varsel)) %>%
  MclustDA(class = pcscores_gal_des_f$taxon_sp, G = 1)
# gal+Des
gal_des_1sp_f <- as_tibble(pcscores_gal_des_f[, -1]) %>%
  dplyr::select(all_of(gal_des_f_varsel)) %>%
  MclustDA(class = rep("gal_des", nrow(pcscores_gal_des_f)), G = 1)
# compare models
sort(c(model2sp = gal_des_2sp_f$bic,  # gal, Des
       model1sp = gal_des_1sp_f$bic)) # gal+Des
# model2sp model1sp 
# 1755.936 1834.842

# Geometric data (males)
data_gal_des_g <- filter(data_gm, taxon_sp == "gallardoi" | taxon_sp == "Deseado")
pca_gal_des_g <- data.frame(droplevels(data_gal_des_g$taxon_sp),
                            gm.prcomp(data_gal_des_g[, -1])$x)
colnames(pca_gal_des_g)[1] <- "taxon_sp"
# variable selection
gal_des_g_varsel <- clustvarsel(pca_gal_des_g[, -1],
                                G = 1:2,
                                search = "greedy",
                                direction = "forward",
                                samp = FALSE)
gal_des_g_varsel <- gal_des_g_varsel$subset
# calculate empirical support for models:
# gal, Des
gal_des_2sp_g <- as_tibble(pca_gal_des_g[, -1]) %>%
  dplyr::select(all_of(gal_des_g_varsel)) %>%
  MclustDA(class = pca_gal_des_g$taxon_sp, G = 1)
# gal+Des
gal_des_1sp_g <- as_tibble(pca_gal_des_g[, -1]) %>%
  dplyr::select(all_of(gal_des_g_varsel)) %>%
  MclustDA(class = rep("gal_des", nrow(pca_gal_des_g)), G = 1)
# compare models
sort(c(model2sp = gal_des_2sp_g$bic,  # gal, Des
       model1sp = gal_des_1sp_g$bic)) # gal+Des
# model2sp model1sp 
# 7425.953 7476.005 

#### vs West cluster (archeforus + scolaroi + zullyae)

# Linear data (males)
data_gal_wes_m <- filter(data_m, taxon_sp == "gallardoi" |
                           taxon_sp == "archeforus" |
                           taxon_sp == "scolaroi" | 
                           taxon_sp == "zullyae")
levels(data_gal_wes_m$taxon_sp)[c(1, 10, 16)] <- "westclust"
pcscores_gal_wes_m <- data.frame(taxon_sp = droplevels(data_gal_wes_m$taxon_sp),
                                 prcomp(data_gal_wes_m[, -1])$x)
# variable selection
gal_wes_m_varsel <- clustvarsel(pcscores_gal_wes_m[, -1],
                                G = 1:2,
                                search = "greedy",
                                direction = "forward",
                                samp = FALSE)
gal_wes_m_varsel <- gal_wes_m_varsel$subset
# calculate empirical support for models:
# gal, Wes
gal_wes_2sp_m <- as_tibble(pcscores_gal_wes_m[, -1]) %>%
  dplyr::select(all_of(gal_wes_m_varsel)) %>%
  MclustDA(class = pcscores_gal_wes_m$taxon_sp, G = 1)
# gal+Wes
gal_wes_1sp_m <- as_tibble(pcscores_gal_wes_m[, -1]) %>%
  dplyr::select(all_of(gal_wes_m_varsel)) %>%
  MclustDA(class = rep("gal_wes", nrow(pcscores_gal_wes_m)), G = 1)
# compare models
sort(c(model2sp = gal_wes_2sp_m$bic,  # gal, Wes
       model1sp = gal_wes_1sp_m$bic)) # gal+Wes
# model2sp model1sp 
# 361.8067 378.8026  

# Linear data (females)
data_gal_wes_f <- filter(data_f, taxon_sp == "gallardoi" |
                           taxon_sp == "archeforus" |
                           taxon_sp == "scolaroi" |
                           taxon_sp == "zullyae")
levels(data_gal_wes_f$taxon_sp)[c(1, 10, 16)] <- "westclust"
pcscores_gal_wes_f <- data.frame(taxon_sp = droplevels(data_gal_wes_f$taxon_sp),
                                 prcomp(data_gal_wes_f[, -1])$x)
# variable selection
gal_wes_f_varsel <- clustvarsel(pcscores_gal_wes_f[, -1],
                                G = 1:2,
                                search = "greedy",
                                direction = "forward",
                                samp = FALSE)
gal_wes_f_varsel <- gal_wes_f_varsel$subset
# calculate empirical support for models:
# gal, Wes
gal_wes_2sp_f <- as_tibble(pcscores_gal_wes_f[, -1]) %>%
  dplyr::select(all_of(gal_wes_f_varsel)) %>%
  MclustDA(class = pcscores_gal_wes_f$taxon_sp, G = 1)
# gal+Wes
gal_wes_1sp_f <- as_tibble(pcscores_gal_wes_f[, -1]) %>%
  dplyr::select(all_of(gal_wes_f_varsel)) %>%
  MclustDA(class = rep("gal_wes", nrow(pcscores_gal_wes_f)), G = 1)
# compare models
sort(c(model2sp = gal_wes_2sp_f$bic,  # gal, Wes
       model1sp = gal_wes_1sp_f$bic)) # gal+Wes
# model2sp model1sp 
# 129.8781 167.2023

# Geometric data (males)
data_gal_wes_g <- filter(data_gm, taxon_sp == "gallardoi" |
                           taxon_sp == "archeforus" |
                           taxon_sp == "scolaroi" | 
                           taxon_sp == "zullyae")
levels(data_gal_wes_g$taxon_sp)[c(1, 10, 16)] <- "westclust"
pca_gal_wes_g <- data.frame(droplevels(data_gal_wes_g$taxon_sp),
                            gm.prcomp(data_gal_wes_g[, -1])$x)
colnames(pca_gal_wes_g)[1] <- "taxon_sp"
# variable selection
gal_wes_g_varsel <- clustvarsel(pca_gal_wes_g[, -1],
                                G = 1:2,
                                search = "greedy",
                                direction = "forward",
                                samp = FALSE)
gal_wes_g_varsel <- gal_wes_g_varsel$subset
# calculate empirical support for models:
# gal, Wes
gal_wes_2sp_g <- as_tibble(pca_gal_wes_g[, -1]) %>%
  dplyr::select(all_of(gal_wes_g_varsel)) %>%
  MclustDA(class = pca_gal_wes_g$taxon_sp, G = 1)
# gal+Wes
gal_wes_1sp_g <- as_tibble(pca_gal_wes_g[, -1]) %>%
  dplyr::select(all_of(gal_wes_g_varsel)) %>%
  MclustDA(class = rep("gal_wes", nrow(pca_gal_wes_g)), G = 1)
# compare models
sort(c(model2sp = gal_wes_2sp_g$bic,  # gal, Wes
       model1sp = gal_wes_1sp_g$bic)) # gal+Wes
# model2sp model1sp 
# 2205.436 2239.547