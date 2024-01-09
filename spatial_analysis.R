# devtools::install_github("danlwarren/ENMTools")
# remotes::install_github("marlonecobos/rangemap")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("treeio")
# devtools::install_github("wxguillo/machuruku")

load("spatial_analysis.RData")
save.image("spatial_analysis.RData")

par(family = "lmsans10")
packages <- c("rangemap", "vcfR", "adegenet", "maptools", "caret", "geoR",
              "varhandle", "raster", "logisticPCA", "fossil", "vegan",
              "machuruku", "ape", "raster", "rworldmap", "terra")
lapply(packages, library, character.only = TRUE) # this loads as many as you put in 'packages'

# Distribution maps ----

ocurr_sp <- read.csv("data/sheets/full_records.csv", header = TRUE)[, -4]
# replace some values that falls out of the continent
ocurr_sp[ocurr_sp == -65.7394444444445] <- -66
ocurr_sp[ocurr_sp == c(-65.8391944444444, -65.8410833333333, -65.8913888888889)] <- -66
ocurr_sp[ocurr_sp == c(-47.7149722222222, -47.7169722222222, -47.7493333333333)] <- -47.7

# check occurrence
par(mar = c(1, 1, 1, 1))
x <- c(-75, -63)
y <- c(-52, -41) # dist. grupo kingii
rangemap_explore(occurrences = ocurr_sp, xlim = x, ylim = y)

# Species ranges from hull polygons

rgeos::set_RGEOS_CheckValidity(2L)
# concave hulls
# do this for each sp:
hull_conc <- rangemap_hull(occurrences = ocurr_sp[ocurr_sp$design_citb_taxon == "archeforus", ],
                           hull_type = "concave",
                           buffer_distance = 20000,
                           extent_of_occurrence = FALSE,
                           area_of_occupancy = FALSE,
                           save_shp = FALSE,
                           name = "data/maps/archeforus_range",
                           overwrite = TRUE)
# full distribution of the kingii group:
hull_conc_full <- rangemap_hull(occurrences = ocurr_sp,
                                hull_type = "concave",
                                buffer_distance = 20000,
                                extent_of_occurrence = FALSE,
                                area_of_occupancy = FALSE,
                                save_shp = FALSE,
                                name = "data/maps/fullrange",
                                overwrite = TRUE)
# plot
rangemap_plot(hull_conc_full, add_occurrences = TRUE)


# Test the effect of eco-geographical predictors on genetic structure ----

#### Correlation between environmental variables and genetic structure using RDA ####
# Based on scripts provided by F. Burbrink (DOI: 10.1111/evo.14141)

# Response variable: Genetic distances

# load genomic data
vcf_full <- read.vcfR("data/genomics/populations.snps.maf0.01.vcf")
gen_full <- vcfR2genind(vcf_full)
ploidy(gen_full) <- 2

# get genetic distances
dist_full <- dist(gen_full@tab)

# specimen mapping to groups and geographic coordinates
data <- read.table("data/sheets/mapping_seqs.csv", header = TRUE, sep = ",")
# specimens LJAMM9474, 11650, and 7352 assigned to "Deseado", according to tess3r

# ensure the matching of rows
ind <- data.frame(indNames(gen_full))
names(ind) <- "ind"
data <- merge(ind, data, by = "ind")
sum(indNames(gen_full) == data$ind) == length(indNames(gen_full)) # everything matches?

# check points on map
data("wrld_simpl")
plot(wrld_simpl,
     xlim = c(min(data$long), max(data$long)),
     ylim = c(min(data$lat), max(data$lat)),
     axes = TRUE, col = "grey95")
points(x = data$long, y = data$lat, col = "black", pch = "+", cex = 0.75)

# localities for everything
locs_all <- cbind(data$long, data$lat)

# Predictor variables:
# geographic distances
# environmental variables (past and current)
# elevation
# physiographic units (vegetation)

# Bioclim data from http://www.paleoclim.org/ (2.5 m resolution)
# you can download the data and put in the bioclim folder within data/

# 1) Pliocene: Marine Isotope Stage M2 (ca. 3.3 Ma), v1.0*,
# 2) Pliocene: mid-Pliocene warm period (3.205 Ma), v1.0*,
# 3) Pleistocene: Marine Isotope Stage M19 (ca. 787 ka), v1.0*,
# 4) Pleistocene: Last Interglacial (ca. 130 ka), v1.0,
# 5) Pleistocene: Last Glacial Maximum (ca. 21 ka), v1.2b**
# 6) Current (1979 – 2013): Anthropocene, v1.2b**

# put all directories in a vector
files <- list.files(pattern = "\\.tif$", recursive = TRUE)
# stack bioclim data
cur <- stack(files[grepl("data/bioclim/cur", files)])
lgm <- stack(files[grepl("data/bioclim/LGM", files)])
lig <- stack(files[grepl("data/bioclim/LIG", files)])
mis19 <- stack(files[grepl("data/bioclim/MIS19", files)])
m2 <- stack(files[grepl("data/bioclim/M2", files)])
mpwp <- stack(files[grepl("data/bioclim/mPWP", files)])

# Altitudinal data of Argentina and Chile

getData("alt", country = "ARG", mask = TRUE, download = TRUE)
getData("alt", country = "CHL", mask = TRUE, download = TRUE)
ele_arg <- raster("ARG_msk_alt.grd")
ele_chi <- raster("CHL1_msk_alt.grd")
ele <- raster::merge(ele_arg, ele_chi)
plot(ele) # check

# pull the bioclim and altitude variables and format each

ele_all <- extract(ele, locs_all)
mpwp_all <- extract(mpwp, locs_all)
lgm_all <- extract(lgm, locs_all)
lig_all <- extract(lig, locs_all)
mis19_all <- extract(mis19, locs_all)
m2_all <- extract(m2, locs_all)
cur_all <- extract(cur, locs_all)

# check there are no NAs in any row
sum(is.na(ele_all)) # elev
which(apply(mpwp_all, 1, function(X) any(is.na(X))), arr.ind = TRUE)
which(apply(lgm_all, 1, function(X) any(is.na(X))), arr.ind = TRUE)
which(apply(lig_all, 1, function(X) any(is.na(X))), arr.ind = TRUE)
which(apply(mis19_all, 1, function(X) any(is.na(X))), arr.ind = TRUE)
which(apply(m2_all, 1, function(X) any(is.na(X))), arr.ind = TRUE)
which(apply(cur_all, 1, function(X) any(is.na(X))), arr.ind = TRUE)

colnames(mpwp_all) <- paste0("mpwp_", colnames(mpwp_all))
colnames(lgm_all) <- paste0("lgm_", colnames(lgm_all))
colnames(lig_all) <- paste0("lig_", colnames(lig_all))
colnames(mis19_all) <- paste0("mis19_", colnames(mis19_all))
colnames(m2_all) <- paste0("m2_", colnames(m2_all))
colnames(cur_all) <- paste0("cur_", colnames(cur_all))

# remove correlated variables (cutoff = 0.9)

lgm_corr <- cor(lgm_all)
lgm_rmv <- findCorrelation(lgm_corr, cutoff = 0.9, verbose = FALSE, names = FALSE, exact = TRUE)
lgm_uncorr <- lgm_all[, -lgm_rmv]

lig_corr <- cor(lig_all)
lig_rmv <- findCorrelation(lig_corr, cutoff = 0.9, verbose = FALSE, names = FALSE, exact = TRUE)
lig_uncorr <- lig_all[, -lig_rmv]

mis19_corr <- cor(mis19_all)
mis19_rmv <- findCorrelation(mis19_corr, cutoff = 0.9, verbose = FALSE, names = FALSE, exact =TRUE)
mis19_uncorr <- mis19_all[, -mis19_rmv]

m2_corr <- cor(m2_all)
m2_rmv <- findCorrelation(m2_corr, cutoff = 0.9, verbose = FALSE, names = FALSE, exact =TRUE)
m2_uncorr <- m2_all[, -m2_rmv]

cur_corr <- cor(cur_all)
cur_rmv <- findCorrelation(cur_corr, cutoff = 0.9, verbose = FALSE, names = FALSE, exact =TRUE)
cur_uncorr <- cur_all[, -cur_rmv]

mpwp_corr <- cor(mpwp_all)
mpwp_rmv <- findCorrelation(mpwp_corr, cutoff = 0.9, verbose = FALSE, names = FALSE, exact =TRUE)
mpwp_uncorr <- mpwp_all[, -mpwp_rmv]

data2 <- data.frame(data,
                    ele_all,
                    lgm_uncorr,
                    lig_uncorr,
                    mis19_uncorr,
                    m2_uncorr,
                    cur_uncorr,
                    mpwp_uncorr)

# add some very small jitter so that overlapping localities are not at the exact same lat/lon
coords_all <- jitterDupCoords(data.frame(data2$long, data2$lat), max = 0.001)
# check distribution of points
plot(wrld_simpl,
     xlim = c(min(coords_all$data2.long), max(coords_all$data2.long)),
     ylim = c(min(coords_all$data2.lat), max(coords_all$data2.lat)),
     axes = TRUE, col = "grey95")
points(x = coords_all$data2.long, y = coords_all$data2.lat, col = "black", pch = "+", cex = 0.75)

# break environmental variables down into PCs
matrix(colnames(data2)) # check names and indexes
lgm_pc <- prcomp(data2[, 7:14], center = TRUE)$x[, 1]
lig_pc <- prcomp(data2[, 15:24], center = TRUE)$x[, 1]
mis19_pc <- prcomp(data2[, 25:32], center = TRUE)$x[, 1]
m2_pc <- prcomp(data2[, 33:38], center = TRUE)$x[, 1]
cur_pc <- prcomp(data2[, 39:45], center = TRUE)$x[, 1]
mpwp_pc <- prcomp(data2[, 46:51], center = TRUE)$x[, 1]

# vegetation units of Argentina
veg_unit <- shapefile("data/maps/UnidadesVegetacionArg_P7F4.shp")
veg_unit <- spTransform(veg_unit, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),
                        use_ob_tran = TRUE)
veg_unit_all <- extract(veg_unit, locs_all)
# dummy code categorical variables
categ_veg <- as.numeric(veg_unit_all$CODIGO)
categ_veg <- to.dummy(categ_veg, "unit")
colnames(categ_veg) <- paste0("veget_", colnames(categ_veg))
# transform these variables into PCA space to break ecoregions down into a single meaningful vector
categ_vegPC <- logisticPCA(categ_veg, k = 2)$PCs[, 1]

# Geographic distances of individuals

geo_dist <- earth.dist(data.frame(coords_all$data2.long, coords_all$data2.lat))
# Principal Coordinates of Neighboring Matrices: Break spatial distances into many eigenvectors
score_dist <- pcnm(geo_dist)$vectors

# use RDA to find which of these spatial variables are significant with respect to genetic distances
# use those in the real RDA analysis
ncol(score_dist)
test_dist <- capscale(dist_full ~ score_dist[, 1] + score_dist[, 2] +
                        score_dist[, 3] + score_dist[, 4] +
                        score_dist[, 5] + score_dist[, 6] +
                        score_dist[, 7] + score_dist[, 8] +
                        score_dist[, 9] + score_dist[, 10] +
                        score_dist[, 11] + score_dist[, 12] +
                        score_dist[, 13] + score_dist[, 14] +
                        score_dist[, 15] + score_dist[, 16] +
                        score_dist[, 17] + score_dist[, 18] +
                        score_dist[, 19] + score_dist[, 20] +
                        score_dist[, 21] + score_dist[, 22] +
                        score_dist[, 23] + score_dist[, 24] +
                        score_dist[, 25] + score_dist[, 26] +
                        score_dist[, 27] + score_dist[, 28] +
                        score_dist[, 29] + score_dist[, 30] +
                        score_dist[, 31] + score_dist[, 32] +
                        score_dist[, 33])

# this will give the significant spatial variables
dist_aov <- anova.cca(test_dist, by = "margin")
pval_dist <- dist_aov$`Pr(>F)`
signif_dist <- which(pval_dist < 0.05)

# rerunning this for significant space variables and environmental variables
# (check there are no NAs)
test_dist2 <- capscale(dist_full ~ categ_vegPC +
                         data2$ele_all +
                         mis19_pc +
                         m2_pc +
                         lgm_pc +
                         lig_pc +
                         cur_pc +
                         mpwp_pc +
                         score_dist[, signif_dist],
                       sqrt.dist = TRUE)
anova.cca(test_dist2, by = "margin")

# Permutation test for capscale under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# Model: capscale(formula = dist_full ~ categ_vegPC + data2$ele_all + mis19_pc + m2_pc + lgm_pc + lig_pc + cur_pc + mpwp_pc + score_dist[, signif_dist], sqrt.dist = TRUE)
#                           Df Variance      F Pr(>F)    
# categ_vegPC                1   0.5558 1.3343  0.012 *  
# data2$ele_all              1   0.5228 1.2550  0.022 *  
# mis19_pc                   1   0.5017 1.2043  0.058 .  
# m2_pc                      1   0.5717 1.3723  0.004 ** 
# lgm_pc                     1   0.4547 1.0914  0.179    
# lig_pc                     1   0.4554 1.0931  0.182    
# cur_pc                     1   0.4804 1.1531  0.088 .  
# mpwp_pc                    1   0.5322 1.2775  0.019 *  
# score_dist[, signif_dist] 13   7.5915 1.4018  0.001 ***
#   Residual                  66  27.4951                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#### Estimation of variable importance with Neural Networks ####

# make a new table suitable for NN
data_nn <- data.frame(categ_vegPC,
                      data2$ele_all,
                      mis19_pc,
                      m2_pc,
                      lgm_pc,
                      mpwp_pc,
                      lig_pc,
                      cur_pc,
                      score_dist[, signif_dist])

# convert genetic distances into PCoA
gendist_pcoa <- dist(gen_full)
gendist_euc <- quasieuclid(gendist_pcoa)
k <- dudi.pco(d = gendist_euc, scannf = FALSE, nf = 10)

# join with the table
data3 <- data.frame(k$l1, data_nn)

# function to perform NN
test_NN_new_pred <- function(data, number, predicted){
  maxs <- apply(data, 2, max) 
  mins <- apply(data, 2, min)
  scaled <- as.data.frame(scale(data, center = mins, scale = maxs - mins))
  data <- scaled
  
  index <- sample(1:nrow(data), round(0.7*nrow(data)))
  train.cv <- data[index, ]
  test.cv <- data[-index, ]
  rem <- setdiff(predicted, number)
  
  train.cv <- train.cv[, -rem]
  test.cv <- test.cv[, -rem]
  # train.cv[, 1] <- as.factor(train.cv[, 1])
  # test.cv[, 1] <- as.factor(test.cv[, 1])
  fit <- train(as.formula(paste(names(train.cv)[1], "~.")),
               data = train.cv,
               method = "nnet",
               maxit = 1000,
               trace = TRUE,
               linout = 1)  
  
  important <- varImp(object = fit)
  imp <- important$importance
  imp2 <- imp[ order(rownames(imp)), , drop = FALSE]
  
  predict <- predict(fit, newdata = test.cv)
  rmse <- sqrt(mean((predict - test.cv[, 1])^2)) 
  
  results <- list()
  results$importance <- imp2
  results$RMSE <- postResample(pred = predict, obs = test.cv[, 1])[1]
  results$r2 <- postResample(pred = predict, obs = test.cv[, 1])[2]
  results$MAE <- postResample(pred = predict, obs = test.cv[, 1])[3]
  new_test.cv <- test.cv
  
  new_test.cv[, 1] <- sample(new_test.cv[, 1])
  
  results$randomRMSE <- postResample(pred = predict, obs = new_test.cv[, 1])[1]
  results$random_r2 <- postResample(pred = predict, obs = new_test.cv[, 1])[2]
  results$randomMAE <- postResample(pred = predict, obs = new_test.cv[, 1])[3]
  
  return(results)
  }

# this will loop the NN analyses for each axis of the PCoA genetic distances
test_nn <- lapply(c(1:10), function(x) test_NN_new_pred(data3, x, c(1:10)))

# this will pull of your outputs from those NN for each axis
imp <- lapply(c(1:10), function(x) test_nn[[x]]$importance)
RMSE <- lapply(c(1:10), function(x) test_nn[[x]]$RMSE)
r2 <- lapply(c(1:10), function(x) test_nn[[x]]$r2)
MAE <- lapply(c(1:10), function(x) test_nn[[x]]$MAE)
randomRMSE <- lapply(c(1:10), function(x) test_nn[[x]]$randomRMSE)
random_r2 <- lapply(c(1:10), function(x) test_nn[[x]]$random_r2)
randomMAE <- lapply(c(1:10), function(x) test_nn[[x]]$randomMAE)

# (plots 6 x 7 inch. ldscape)
# plot your accuracy on each PCoA genetic distance axis against random accuracy
plot(c(1:10), r2, col = "blue", pch = 19, cex = 2,
     xlab = "PCoA axis", ylab = "accuracy", cex.lab = 1.5)
points(c(1:10), random_r2, col = "red", pch = 19, cex = 2)

# this function will plot variable importance for each axis
plot_bar <- function(data, x){
  j <- data[[x]]
  s <- j$Overall
  names(s) <- row.names(j)
  barplot(rev(sort(s)), cex.names = 0.8, width = 0.5, las = 2)
  return(rev(sort(s)))
  }

# variable importance of the first three axes
plot_bar(imp, 1)
plot_bar(imp, 2)
plot_bar(imp, 3)


#### Spatial and ecological divergence btw kingii and Deseado clade ####
# since these two populations diverged ~150Kya,
# we will only use current spatial/climate data,
# LGM, and LIG climate data 

vcf_kindes <- read.vcfR("data/genomics/populations.snps.maf0.01.kin_des.vcf")
gen_kindes <- vcfR2genind(vcf_kindes)
ploidy(gen_kindes) <- 2
# get genetic distances
dist_kindes <- dist(gen_kindes@tab)

data_kindes <- data[data$ind %in% indNames(gen_kindes), ]

locs_kindes <- data.frame(data_kindes$long, data_kindes$lat)
ele_kindes <- extract(ele, locs_kindes)
lgm_kindes <- extract(lgm, locs_kindes)
lig_kindes <- extract(lig, locs_kindes)
cur_kindes <- extract(cur, locs_kindes)

veg_unit_kingdes <- extract(veg_unit, locs_kindes)
# dummy code categorical variables
categ_veg_kindes <- as.numeric(veg_unit_kingdes$CODIGO)
categ_veg_kindes <- to.dummy(categ_veg_kindes, "unit")
colnames(categ_veg_kindes) <- paste0("veget_", colnames(categ_veg_kindes))
# transform these variables into PCA space to break ecoregions down into a single meaningful vector
categ_vegPC_kindes <- logisticPCA(categ_veg_kindes, k = 2)$PCs[, 1]

# check there are no NAs in any row
sum(is.na(ele_kindes)) # elev
which(apply(lgm_kindes, 1, function(X) any(is.na(X))), arr.ind = TRUE)
which(apply(lig_kindes, 1, function(X) any(is.na(X))), arr.ind = TRUE)
which(apply(cur_kindes, 1, function(X) any(is.na(X))), arr.ind = TRUE)

colnames(lgm_kindes) <- paste0("lgm_", colnames(lgm_kindes))
colnames(lig_kindes) <- paste0("lig_", colnames(lig_kindes))
colnames(cur_kindes) <- paste0("cur_", colnames(cur_kindes))

# remove correlated variables (cutoff = 0.9)

lgm_corr_kindes <- cor(lgm_kindes)
lgm_rmv_kindes <- findCorrelation(lgm_corr_kindes, cutoff = 0.9, verbose = FALSE, names = FALSE, exact = TRUE)
lgm_uncorr_kindes <- lgm_kindes[, -lgm_rmv_kindes]

lig_corr_kindes <- cor(lig_kindes)
lig_rmv_kindes <- findCorrelation(lig_corr_kindes, cutoff = 0.9, verbose = FALSE, names = FALSE, exact = TRUE)
lig_uncorr_kindes <- lig_kindes[, -lig_rmv_kindes]

cur_corr_kindes <- cor(cur_kindes)
cur_rmv_kindes <- findCorrelation(cur_corr_kindes, cutoff = 0.9, verbose = FALSE, names = FALSE, exact =TRUE)
cur_uncorr_kindes <- cur_kindes[, -cur_rmv_kindes]

data2_kd <- data.frame(data_kindes, ele_kindes, lgm_uncorr_kindes, lig_uncorr_kindes, cur_uncorr_kindes)

# add some very small jitter so that overlapping localities are not at the exact same lat/lon
coords_kindes <- jitterDupCoords(data.frame(data2_kd$long, data2_kd$lat), max = 0.001)
# check distribution of points
plot(wrld_simpl,
     xlim = c(min(coords_kindes$data2_kd.long), max(coords_kindes$data2_kd.long)),
     ylim = c(min(coords_all$data2.lat), max(coords_all$data2.lat)),
     axes = TRUE, col = "grey95")
points(x = coords_kindes$data2_kd.long,
       y = coords_kindes$data2_kd.lat,
       col = "black", pch = "+", cex = 0.75)

# break environmental variables down into PCs
matrix(colnames(data2_kd)) # check names and indexes
lgm_pc_kd <- prcomp(data2_kd[, 7:15], center = TRUE)$x[, 1]
lig_pc_kd <- prcomp(data2_kd[, 16:23], center = TRUE)$x[, 1]
cur_pc_kd <- prcomp(data2_kd[, 24:31], center = TRUE)$x[, 1]

# Geographic distances of individuals
geo_dist_kd <- earth.dist(data.frame(coords_kindes$data2_kd.long, coords_kindes$data2_kd.lat))
# Principal Coordinates of Neighboring Matrices: Break spatial distances into many eigenvectors
score_dist_kd <- pcnm(geo_dist_kd)$vectors

# use RDA to find which of these spatial variables are significant with respect to genetic distances
# use those in the real RDA analysis
ncol(score_dist_kd)
test_dist_kd <- capscale(dist_kindes ~ score_dist_kd[, 1] +
                           score_dist_kd[, 2] + score_dist_kd[, 3] +
                           score_dist_kd[, 4] + score_dist_kd[, 5] +
                           score_dist_kd[, 6] + score_dist_kd[, 7] +
                           score_dist_kd[, 8] + score_dist_kd[, 9] +
                           score_dist_kd[, 10] + score_dist_kd[, 11] +
                           score_dist_kd[, 12] + score_dist_kd[, 13] +
                           score_dist_kd[, 14] + score_dist_kd[, 15])

# this will give the significant spatial variables
dist_aov_kd <- anova.cca(test_dist_kd, by = "margin")
pval_dist_kd <- dist_aov_kd$`Pr(>F)`
signif_dist_kd <- which(pval_dist_kd < 0.05)
# no significant spatial variables

# here I rerun considering all the variables,
# but given that there are no significant spatial variables, I do not subset as in the previous analysis  
test_dist2_kd <- capscale(dist_kindes ~ categ_vegPC_kindes +
                            data2_kd$ele_kindes +
                            lgm_pc_kd +
                            lig_pc_kd +
                            cur_pc_kd +
                            score_dist_kd[, signif_dist_kd],
                          sqrt.dist = TRUE)
anova.cca(test_dist2_kd, by = "margin")

# Permutation test for capscale under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# Model: capscale(formula = dist_kindes ~ categ_vegPC_kindes + data2_kd$ele_kindes + lgm_pc_kd + lig_pc_kd + cur_pc_kd + score_dist_kd[, signif_dist_kd], sqrt.dist = TRUE)
#                                 Df Variance      F Pr(>F)  
# categ_vegPC_kindes               1    1.101 0.9434  0.913  
# data2_kd$ele_kindes              1    1.130 0.9676  0.728  
# lgm_pc_kd                        1    1.136 0.9729  0.636  
# lig_pc_kd                        1    1.135 0.9725  0.657  
# cur_pc_kd                        1    1.133 0.9702  0.670  
# score_dist_kd[, signif_dist_kd]  2    2.553 1.0934  0.044 *
# Residual                        36   42.031                
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Estimation of variable importance

data_kd_nn <- data.frame(categ_vegPC_kindes,
                         data2_kd$ele_kindes,
                         lgm_pc_kd, lig_pc_kd,
                         cur_pc_kd,
                         pcnm(geo_dist_kd)$vectors)
# convert genetic distances into PCoA
gendist_pcoa_kd <- dist(gen_kindes)
gendist_euc_kd <- quasieuclid(gendist_pcoa_kd)
k_kd <- dudi.pco(d = gendist_euc_kd, scannf = FALSE, nf = 10)
data3_kd <- data.frame(k_kd$l1, data_kd_nn)
# NN analyses
test_nn_kd <- lapply(c(1:10), function(x) test_NN_new_pred(data3_kd, x, c(1:10)))
# pull of your outputs from those NN for each axis
imp_kd <- lapply(c(1:10), function(x) test_nn_kd[[x]]$importance)
RMSE_kd <- lapply(c(1:10), function(x) test_nn_kd[[x]]$RMSE)
r2_kd <- lapply(c(1:10), function(x) test_nn_kd[[x]]$r2)
MAE_kd <- lapply(c(1:10), function(x) test_nn_kd[[x]]$MAE)
randomRMSE_kd <- lapply(c(1:10), function(x) test_nn_kd[[x]]$randomRMSE)
random_r2_kd <- lapply(c(1:10), function(x) test_nn_kd[[x]]$random_r2)
randomMAE_kd <- lapply(c(1:10), function(x) test_nn_kd[[x]]$randomMAE)
# plot your accuracy on each PCoA genetic distance axis against random accuracy
plot(c(1:10), r2_kd, col = "blue", pch = 19, cex = 2,
     xlab = "PCoA axis", ylab = "accuracy", cex.lab = 1.5)
points(c(1:10), random_r2_kd, col = "red", pch = 19, cex = 2)

# variable importance of the first three axes
plot_bar(imp_kd, 7)
plot_bar(imp_kd, 8)
plot_bar(imp_kd, 3)


## dada una elevada influencia de variables espaciales en la divers. genetica, probar estadisticamente IBD
# mantel.randtest()



# Phylogenetic niche modelling ----

# load tree and occurrence data
tree <- read.nexus("data/phylo/mcc_median.tre")
occ	<- read.delim("data/sheets/full_records.csv", header = TRUE, sep = ",")
occ <- occ[, -c(1, 5)]; colnames(occ)[1] <- "species"
# account for spatial autocorrelation by rarefying occurrence data
occ_sub <- machu.occ.rarefy(in.pts = occ, rarefy.dist = 5, verbose = FALSE, plot = FALSE)
# visualize points
plot(getMap(resolution = "low"), xlim = c(-77, -57),
     ylim = c(-57, -40), axes = TRUE, col = "grey95")
points(x = occ_sub$long, y = occ_sub$lat, col = "black", pch = "+", cex = 0.75)
# check each group's data points
summary(factor(occ_sub$species))

# load climate data
cur <- rast(list.files("data/biclim/cur", full.names = TRUE, pattern = "\\.tif$"))
lgm <- rast(list.files("data/biclim/LGM", full.names = TRUE, pattern = "\\.tif$"))
lig <- rast(list.files("data/biclim/LIG", full.names = TRUE, pattern = "\\.tif$"))
mis19 <- rast(list.files("data/biclim/MIS19", full.names = TRUE, pattern = "\\.tif$"))
mpwp <- rast(list.files("data/biclim/mPWP", full.names = TRUE, pattern = "\\.tif$"))
m2 <- rast(list.files("data/biclim/M2", full.names = TRUE, pattern = "\\.tif$"))

# delimit extent of layers
limits <- c(-77, -57, -57, -40)
ext <- as(extent(limits), "SpatialPolygons")
crs(ext) <- "+proj=longlat +datum=WGS84 +no_defs"
cur_dl <- crop(cur, ext)
lgm_dl <- crop(lgm, ext)
lig_dl <- crop(lig, ext)
mis19_dl <- crop(mis19, ext)
m2_dl <- crop(m2, ext)
mpwp_dl <- crop(mpwp, ext)

# Select variables with greatest contributions to species’ niches
# It's important to make sure that all time periods have the same set of climate variables
cur_var <- machu.top.env(occ_sub, cur_dl, method = "nvars", nvars.save = 4, verbose = TRUE)
cur_var <- c("bio_15", "bio_12", "bio_18", "bio_8")
# the following variables are not included in paleoclimate layers mis19, m2, and mpwp:
# (you can omit this line if the selected variables does not include any of the following)
# cur_var <- setdiff(cur_var, c("bio_2", "bio_3", "bio_5", "bio_6", "bio_7"))
cur_sub_dl <- cur_dl[[cur_var]]
lgm_sub_dl <- lgm_dl[[cur_var]]
lig_sub_dl <- lig_dl[[cur_var]]
mis19_sub_dl <- mis19_dl[[cur_var]]
m2_sub_dl <- m2_dl[[cur_var]]
mpwp_sub_dl <- mpwp_dl[[cur_var]]

#### 1. Estimating tip response curves ####

resp <- machu.1.tip.resp(occ_sub, cur_sub_dl, verbose = TRUE, plot = "f", plot.points = FALSE)
resp <- resp[[2]]
# visualizing climate response curves
machu.respplot(resp, plot = "t")

#### 2. Estimating ancestral niches at time slices ####

# plot time slices
# machu.treeplot(tree, timelaboffset = -.1, timeslice = c(0, 0.021, 0.13, 0.787, 3.205, 3.3))

# Estimate ancestral niches at timeslices
# accounting for uncertainty in ancestral character estimation
ace_ts <- machu.2.ace(resp, tree, timeslice = c(0.021, 0.787, 3.205), unc = TRUE)
# visualize uncertainty
# machu.respplot(ace_ts[[1]], clim = "bio_10")
# also estimate niches for all nodes
ace_all <- machu.2.ace(resp, tree, unc = TRUE)

# visualize how these taxa (current and ancestral) have evolved in climate space
# machu.respplot(ace_all$tips_and_nodes, clim = "bio_18")

#### 3. Projecting ancestral models into paleoclimatic data ####

# combine paleoclimates in a list
# ensure the datasets are listed in the same order as their corresponding timeslices in the machu.2.ace output
clim <- list(lgm_sub_dl, lig_sub_dl, mis19_sub_dl, mpwp_sub_dl, m2_sub_dl)
model_ts <- machu.3.anc.niche(ace_ts, clim, verbose = TRUE)
# build present-day models for only extant taxa
model_cur <- machu.3.anc.niche(ace_all, cur_sub_dl, taxa = 1:length(tip.label(sp_tree)), verbose = TRUE)
# plot
machu.plotmap(model_ts, plot = "t", axes = FALSE, title.cex = 0.85, col = viridis::magma(100))

# In case you want to export the rasters
dir.create("ENMS_rasters")
machu.3.anc.niche(ace_ts, clim, output.folder = "ENMs_rasters", verbose = TRUE) %>% invisible

#### Plots ####
map <- shapefile("data/maps/ne_50m_admin_0_countries.shp")
map_crop <- crop(x = map, y = limits)

dir.create("ENMs_imgs")
svg("ENMs_imgs/m2_kinClade", width = 5, height = 7)
par(family = "Latin Modern Sans")
plot(model_ts$timeslice_3.3$`Node1-Node2`, col = viridis::magma(100), axes = FALSE, legend = FALSE)
lines(map_crop, col = "white", lty = "dashed", lwd = 1.417)
# lines(shapefile("../GIS/rangemap/somuncurae_range.shp"), col = "white", lwd = 3.78)
# lines(shapefile("../GIS/rangemap/uptoni_range.shp"), col = "white", lwd = 3.78)
# lines(shapefile("../GIS/rangemap/zullyae_range.shp"), col = "white", lwd = 3.78)
# lines(shapefile("../GIS/rangemap/sarmientoi_range.shp"), col = "white", lwd = 3.78)
dev.off()
