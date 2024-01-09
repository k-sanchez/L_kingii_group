# remotes::install_github("eriqande/TESS3_encho_sen", force = TRUE) # special version of tess3r
# remotes::install_github("eriqande/genoscapeRtools") # plot genoscapes
# devtools::install_github("royfrancis/pophelper") # manipulate tess3r output
# devtools::install_github("bcm-uga/LEA")
# devtools::install_github("dipetkov/eems/plotting/rEEMSplots")
# devtools::install_github("dosreislab/bppr")
packages <- c("adegenet", "dartR", "vcfR", "RColorBrewer", "tess3r",
              "pophelper", "sf", "raster", "genoscapeRtools", "ggspatial",
              "phangorn", "LEA", "bppr", "coda", "pegas", "delimitR", "ape",
              "phytools", "TreeTools", "ggplot2")
lapply(packages, library, character.only = TRUE) # this loads as many as you put in "packages"
par(family = "lmsans10") # I use this font in the plots


# I recommend to save the objects in a .Rdata because some of them might take some time to be outputted
load("phylogeography.RData")
save.image("phylogeography.RData")



# Spatial genomic structure ----

vcf_full_file <- "data/genomics/populations.snps.maf0.01.vcf"
vcf_full <- read.vcfR(vcf_full_file) # 11491 SNPs
gi_full <- vcfR2genind(vcf_full)
ploidy(gi_full) <- 2

#### tess3r ####

dir.create("tess3r")
setwd("tess3r")
gl_full <- gi2gl(gi_full)
coords <- as.data.frame(read.csv("../data/sheets/mapping_seqs.csv"))
tess3_full <- tess3(X = as.data.frame(gl_full),
                    coord = as.matrix(coords[, 4:5]),
                    K = 1:15, method = "projected.ls",
                    ploidy = 2, rep = 100, openMP.core.num = 4)
plot(tess3_full, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score",
     main = "Cross-validation criterion") # 5x5 pdf

#### LEA ####

dir.create("../LEA")
setwd("../LEA")
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R") # conversion between formats
vcf2lfmm(vcf_full_file)
# be aware that the output of vcf2lfmm() will be located in the same folder as the input
# if you have .ped data:
# ped2geno("populations.snps.maf0.05.ped")

# Inference of individual admixture coefficients using SNMF
# Main options:
# K = number of ancestral populations
# entropy = TRUE: computes the cross-entropy criterion,
# CPU = the number of CPUs
project <- NULL
project <- snmf("../data/genomics/populations.snps.maf0.01.geno",
                K = 1:15, ploidy = 2, entropy = TRUE, repetitions = 100,
                project = "new", alpha = 10000, CPU = 4)
# to continue a previously interrupted run -> project = "continue"
# plot cross-entropy criterion for all runs in the SNMF project
project <- load.snmfProject("populations.snps.maf0.01.snmfProject")
plot(project, col = "blue", pch = 19, cex = 1.2,
     cex.lab = 1.2, cex.axis = 1.2, main = "Cross-entropy criterion") # pdf 5x5
# select the best run for K
best <- which.min(cross.entropy(project, K = 8))

# barplot
colors <- RColorBrewer::brewer.pal(8, "Set2") # qualitative
bp <- LEA::barchart(project, K = 8, run = best, border = NA, space = 0,
                    col = colors,
                    xlab = "Individuals",
                    ylab = "Ancestry proportions",
                    main = "Ancestry matrix")
axis(1, at = 1:length(bp$order), labels = coords$ind[bp$order],
     las = 1, cex.axis = .4) # plot pdf 4 x 20 lscape

#### snapclust ####

k_aic <- snapclust.choose.k(15, gi_full, IC = AIC)
plot(k_aic, type = "b", cex = 2, xlab = "k", ylab = "AIC")
points(which.min(k_aic), min(k_aic), col = "blue", pch = 20, cex = 2)
abline(v = 8, lty = 2, col = "red")

res6 <- snapclust(gi_full, k = 6)
res7 <- snapclust(gi_full, k = 7)
res8 <- snapclust(gi_full, k = 8)
compoplot(res8, show.lab = TRUE, cex.label = 0.1) # plot pdf 4 x 20 lscape

#### Genomic landscape with optimal K ####

setwd("../tess3r")
# retrieve tess3 Q matrix for K clusters
q_matrix <- qmatrix(tess3_full, K = 8, rep = "best")
write.csv(q_matrix, "ancestry_inds_K8.csv", row.names = coords$ind,
          quote = FALSE)

# STRUCTURE-like barplot for the Q-matrix
# ordered according to the ML tree
tree <- phytools::read.newick("../data/phylo/out_all.treefile")
qlist_tess3r <- readQTess3(tess3_full)
for(i in seq_along(qlist_tess3r)){
  rownames(qlist_tess3r[[i]]) <- coords$ind
  }
qlist_lea8 <- data.frame(Q(project, K = 8,
                           run = which.min(cross.entropy(project, K = 8))))
rownames(qlist_lea) <- coords$ind
qlist_snapclust <- data.frame(res8$proba)
qlist <- list(tess3r = qlist_tess3r$sample8,
              lea = qlist_lea,
              snapclust = qlist_snapclust)
qlist_reorder <- lapply(qlist, function(x) {
  x[match(tree$tip.label[-(1:2)], rownames(coords$ind)), ]
  }
                        )
p <- plotQ(as.qlist(qlist_reorder), imgoutput = "join",
           returnplot = TRUE, exportplot = FALSE, basesize = 5,
           clustercol = colors, useindlab = TRUE, showindlab = TRUE,
           font = "lmsans10")
p # lscape 20x10

# genoscape
full_range <- st_read("../data/maps/fullrange_crop.shp")
ggplot(full_range) + geom_sf() + theme_bw() # check geographic range
# interpolate Q-values by Kriging
# (agrego un "ind." del cluster South en una localidad cerca de la desembocadura del Chalía-Chico)
q_matrix_ed <- as.data.frame(q_matrix)
q_matrix_ed[nrow(q_matrix_ed) + 1, ] <- q_matrix_ed[86, ]
coords_ed <- coords[, 4:5]
coords_ed[nrow(coords_ed) + 1, ] <- c(-69.14756, -50.27572)
names(colors) <- paste("pop", 1:8, sep = "")
genoscape_brick <- tess3r::tess3Q_map_rasters(
  x = q_matrix_ed, 
  coord = coords_ed,
  col.palette = colors,
  map.polygon = full_range,
  window = extent(full_range)[1:4],
  resolution = c(600, 600), # if you want more cells in your raster, set higher
  method = "map.max", 
  interpol = FieldsKrigModel(10),  
  main = "Ancestry coefficients",
  xlab = "Longitude", 
  ylab = "Latitude", 
  cex = .4)

# add names of the clusters back onto this raster brick
names(genoscape_brick) <- c(paste("pop", 1:ncol(q_matrix), sep = ""))

# scaling and cleaning the genoscape_brick
genoscape_rgba <- qprob_rando_raster(
  TRB = genoscape_brick,
  cols = colors, 
  alpha_scale = 4, 
  abs_thresh = 0.0, 
  alpha_exp = 1.55, 
  alpha_chop_max = 230)

# add projection info
crs(genoscape_rgba) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# plot genoscape
ggplot() + ggspatial::layer_spatial(genoscape_rgba) + theme_void() + coord_sf()

# country boundaries and coastlines
boundaries <- st_read("../data/maps/ne_50m_admin_0_countries.shp")
range <- c(xmin = extent(full_range)[1] - 2.5,
           xmax = extent(full_range)[2] + 2.5,
           ymin = extent(full_range)[3] - 1,
           ymax = extent(full_range)[4] + 1) # set range
coast_cropped <- st_crop(boundaries, range) # crop
ggplot(coast_cropped) + geom_sf() + coord_sf() # check cropped part

# plot under default lat-long projection
mapg <- ggplot() + geom_sf(data = coast_cropped) +
layer_spatial(genoscape_rgba) +
theme_void(base_family = "LM Sans 10")
mapg +
coord_sf() +
geom_spatial_point(data = coords_ed, mapping = aes(x = long, y = lat)) # PDF landscape 10 x 8



# SplitsTree network ----
# with ancestry proportions as pies

nnet <- read.nexus.networx("../data/phylo/NNet_uncorrP.nex") # load network
q_matrixDF <- as.data.frame(q_matrix); rownames(q_matrixDF) <- coords$ind
ancestry <- q_matrixDF[nnet$tip.label, ] # sort ancestries in order of tip labels

plot(nnet, "2D", cex = 0.5)
tiplabels(pie = ancestry, piecol = colors, cex = 0.8) # add pies
# export PDF: Landscape 10 x 10


# Effective migration surfaces ----

dir.create("../EEMS")
setwd("../EEMS")
# export PDF
dir.create("plots")
rEEMSplots::eems.plots(mcmcpath = c("nDemes500_chain1/",
                                    "nDemes500_chain2/",
                                    "nDemes500_chain3/"),
                       longlat = TRUE,
                       projection.in = "+proj=longlat +datum=WGS84",
                       projection.out = "+proj=merc +datum=WGS84",
                       add.grid = TRUE,
                       add.map = TRUE,
                       add.demes = TRUE,
                       col.demes = "black",
                       col.map = "gray",
                       lwd.map = 1,
                       out.png = FALSE,
                       plotpath = "plots")




# Concordance Factors table for network estimation ----
# here I use the SNPs2CF() function

# source the function:
# https://github.com/melisaolave/SNPs2CF
source("/directory/of/package/SNPs2CF/functions_v1.6.R")

dir.create("../phylonetworks")
setwd("../phylonetworks")
# Opcional: Convert phased VCF to phy
inputVCF <- "../data/genomics/populations.snps.maf0.01.net.vcf"
vcf <- read.vcf(inputVCF, to = 100000)
vcf
vcf2phylip(vcf.name = inputVCF, total.SNPs = dim(vcf)[2],
           output.name = "populations.snps.maf0.01.phy")
inputPHY <- "populations.snps.maf0.01.phy" # PUT THE .phy MATRIX IN YOUR WORKING DIRECTORY

# Generate CF table
# for Network estimation: n.quartets = 100, bootstrap = TRUE
# for goodness of fit tests: n.quartets = 1, bootstrap = FALSE
output <- SNPs2CF(seqMatrix = inputPHY,
                  ImapName = "../data/sheets/popmap.txt",
                  between.sp.only = TRUE,
                  max.SNPs = NULL,
                  bootstrap = TRUE,
                  outputName = "CF_btw_sp.csv",
                  save.progress = FALSE,
                  n.quartets = 100,
                  # rm.outgroup = TRUE,
                  # outgroupSp = "lineomaculatus",
                  indels.as.fifth.state = FALSE,
                  boots.rep = 100,
                  starting.sp.quartet = 1,
                  max.quartets = 100000,
                  cores = 4)



# delimitR ----

setwd("delimitR")
# to install, set your working directory to the directory above the delimitR directory
# setwd("/media/kevin/KEVIN_HDD/software/paquetes_R/")
# devtools::install("delimitR", dependencies = TRUE)
# .libPaths("~/software/R_packages") # I use this command in our lab computer

# move to the corresponding folder and perform the analysis
setwd("mig_att_som/") # test migration between attenboroughi and somuncurae_uptoni
setwd("mig_bag_sou/") # test migration between baguali and the "South" group
setwd("mig_cha_wes/") # test migration between chacabucoense and the "West" group
setwd("mig_kin_des/") # test migration between kingii and the "Deseado" clade

# Sys.setenv(PATH = paste("/directory/of/package/delimitR/", Sys.getenv("PATH"), sep = ":")) # trick to call python2
Sys.setenv(PATH = paste("~/software/miniconda3/envs/py2_env/bin/",
                        Sys.getenv("PATH"), sep = ":")) # the same in out lab computer
system("python --version") # check python 2

# SFSs generated with easySFS:
# https://github.com/isaacovercast/easySFS
observedSFS <- "input/fastsimcoal2/populations_MSFS" # DO NOT USE THE EXTENSION
traitsfile <- "input/popmap.txt" # numbers AFTER downsampling

#### Building models ####

# An important consideration is population order. Since we need to ensure that our observed SFS is built in the same way as our simulated SFS,
# we must ensure that population order is consistent between our observed SFS and the priors we set up for the analysis.
# The following rules must be followed:
#   1. Populations must be named with numbers 0 to n-1, where n is the number of populations
#   2. For the purposes of building the observed SFS, populations need to be ordered from population 0 to population n-1
# ?setup_fsc2

# Guide tree
# you can also supply a list of guide trees -> list("tree1", "tree2")
obstree <- "(0,1);" # (pop0,pop1)

# Information

# dimensionality -> MAX number of tested species
obsspecies <- 2

# sample sizes and number of SNPs:
# numbers of haploid individuals AFTER downsampling, from sp 0 to sp n-1
# att, som: 6, 6
# bag, sou: 16, 12
# cha, wes: 6, 16
# kin, des: 24, 48
obssampsize <- c(24, 48)
# number of linkage blocks to simulate; for unlinked SNPs, this is the number of SNPs used to BUILD observed SFS:
# att, som: 10262
# bag, sou: 10428
# cha, wes: 15101
# kin, des: 18497
obssnps <- 18497

# Migration

# migration matrix (symmetric)
migmatrix <- matrix(c(FALSE, TRUE,
                      TRUE,  FALSE),
                    nrow = 2, ncol = 2, byrow = TRUE)

# IMscondary contact and divergence with gene flow
IMsccontact <- TRUE    # gene flow begins in the present and ends half-way to the first coalescent event
divwgeneflow <- TRUE # only between sister species; gene flow will begin half-way between time 0 and the coalescent event involving the two sister species, and will end when the two sister species coalesce

# maximum number of migration events in a model,
# must correspond with matrix
maxedges <- 1

# pop sizes in terms of number of haploid individuals (Ne * 2)
# uniform dist. (c(min, max))
obspopsizeprior <- list(2*c(1e4, 1e6),
                        2*c(1e4, 1e6)) # broad priors

# divergence times in number of generations
# divide number of years since divergence by two (1 generation)
# uniform dist. (c(min, max))
# should be provided in order of coalescent intervalthe default model set is inadequate user might generate the fsc26 input files and place them in the working directory

obsdivtimeprior <- list(c(1e4, 1e5)) # broad priors

# coalescence times rules
# coal_rules <- c("Tdiv2$>Tdiv1$")

# migration rates m (one prior for all migration rates) -> ojo! no es lo que arroja 3s (M  = N*m)
# uniform dist. (c(min, max))
obsmigprior <- c(1e-6, 1e-4) # broad priors (= 0.01-100 migrants/gen)

# output prefix:
# name for fsc2 input files
obsprefix <- "div_SC_IM" # divergence, sec. contact, and IM
# should be unique for all guide tree + prior combinations in a folder, or files will be overwritten

#### Setting up the default model set ####
# generate fsc26 input files
# if the default model set is inadequate user might generate the fsc26 input files and place them in the working directory

setup_fsc2(tree = obstree,
           nspec = obsspecies,
           samplesizes = obssampsize,
           nsnps = obssnps,
           prefix = obsprefix,
           migmatrix = migmatrix,
           popsizeprior = obspopsizeprior,
           divtimeprior = obsdivtimeprior,
           migrateprior = obsmigprior,
           secondarycontact = IMsccontact,
           divwgeneflow = divwgeneflow,
           maxmigrations = maxedges)

# If you want to remove a specific model, e.g., a one-species model, delete the corresponding .tpl and .est files and move the remaining model outputs to positions 1, 2, 3, and so on

#### Building a prior ####

# Simulating data

fastsimcoalsims(prefix = obsprefix,
                pathtofsc = "~/fsc26", # the path to the fastsimcoal2 executable
                nreps = 10000) # 10K

# Binning strategy to further summarize the mSFS

nclasses <- min(obssampsize)
# ideally, this number should not be greater than the sample size of the population with the fewest samples
# large values lead to a more complete summary of the data, but also lead to a more sparsely sampled SFS and increased computation times
fullPrior <- makeprior(prefix = obsprefix,
                       nspec = obsspecies,
                       nclasses = nclasses,
                       mydir = getwd(),
                       traitsfile = traitsfile,
                       threshold = 100,
                       thefolder = "prior",
                       ncores = 18)

clean_working(prefix = obsprefix) # remove unnecessary files/folders
reducedPrior <- Prior_reduced(fullPrior) # rm rows with zero variance

# Reading and preparing observed data (binning strategy to observed data)

myobserved <- prepobserved(
  observedSFS,
  fullPrior,
  reducedPrior,
  nclasses,
  obsspecies,
  traitsfile = traitsfile,
  threshold = 100)

#### Building a RF classifier ####
# Predictor variables -> bins of the bSFS
# Response variable -> model used to simulate the data

myRF <- RF_build_abcrf(reducedPrior, fullPrior, 10000)
myRF
# plot(myRF1, training = reducedPrior1)

#### Selecting a model ####

# Applying the RF classifier to the observed data

prediction <- RF_predict_abcrf(myRF, myobserved, reducedPrior, fullPrior, 10000)
prediction
# save.image(file = "delimitR.RData")

# write results to file
write.csv(as.matrix(prediction), file = "prediction.csv")

# write out model info to screen
files <- list.files(pattern = "*.tpl")

for (i in files){
  cat(i, readLines(i), "\n", "\n")
}

# histogram of number of votes vs model (I didn't use this)
# df1 <- data.frame(model = c(1:4), Nvotes = c(0, 152, 207, 141))
# df2 <- data.frame(model = c(1:13), Nvotes = c(0, 138, 33, 111, 13, 129, 0, 1, 28, 21, 0, 1, 25))
# df3 <- data.frame(model = c(1:13), Nvotes = c(0, 0, 10, 1, 17, 5, 52, 48, 77, 5, 149, 85, 11))
# ggplot(df1) + geom_bar(aes(x = model, y = Nvotes), stat = "identity") +
#   theme_bw(base_size = 15, base_family = "LM Sans 10") + 
#   labs(x = "Model", y = "N votes") + 
#   scale_x_continuous("Model", labels = as.character(df3$model), breaks = df3$model)



# Savage-Dickey density ratio ----
# Bayes Factor approximation to assess the support for introgression models
# more details are provided in the Supplementary Material (Methods S2)
setwd("bpp_A00_MSci/step2")

# plot the backbone tree to check the names of the branches involved in introgressions:
plot.phylo(read.tree(text = "((att,nor)A,((bag,sou)B,((wes,cha)C,(des,kin)D)F)G)R;"), show.node.label = TRUE)

# move to the corresponding folder and perform the analysis
setwd("introg1_bid_desD-kinD") # this is the only case of bidirectional introgression
setwd("introg2_wesC-chaC")
setwd("introg3_souB-bagB")
setwd("introg4_CF-attA")
setwd("introg5_chaC-attA")


# load posterior distribution of parameters
table <- read.table(dir(pattern = "mcmc")[1], header = TRUE)
head(table)

# explore mcmc
mcmc.summary(table)
plot(density(table$phi_Y4..X4), col = "blue")

# compare densities of hybrid node time and divergence time of the nearest node (parental of the recipient of the introgression)
plot(density(table$tau_4X), col = "blue", ylim = c(0, 5000), xlim = c(0, 0.005))
lines(density(table$tau_3D), col = "red")
data.frame(mean = round(mcmc.summary(table)$means, 5))

0.0005 / (sum(table$phi_Y2..X2 < 0.0005) / nrow(table))
0.005 / (sum(table$phi_Y2..X2 < 0.005) / nrow(table))
0.05 / (sum(table$phi_Y2..X2 < 0.05) / nrow(table))



# Network calibration ----
# Time calibration of MSCi model

setwd("bpp_A00_MSci/step3/") # folder that includes the MCMC parameters of the full MSCi model  
mcmc_msci <- read.table("mcmc_comb.txt", header = TRUE)

# define function to replace the nth character or string
str_replace_n <- function(x, pattern, replace, n){
  g <- gregexpr(pattern, x)[[1]][n]
  output <- paste0(substr(x, 1, g - 1), replace, substr(x, g + 1, nchar(x)))
  output
}

# Calibrate using a prior on molecular substitution rate:
# this converts tau into absolute divergence times
# Here I use the mean neutral substitution rate estimated for lizards (substitutions/site/Myears)
# taken from Perry et al. 2018 (DOI: 10.1093/gbe/evy157):
umeanmy <- mean(c(0.00056, 0.00058, 0.00069)) # glass lizard, bearded dragon, and anole
# convert to per-generation mutation rate (1 gen ≈ 2 years)
umean <- umeanmy*2/1e6
usd <- 1e-10

# generation time prior
gmean <- 2
gsd <- 0.1

# Explore gamma densities
# calculate the shape (alpha) and scale (beta) parameters
shape <- (umean / usd) ^ 2
scale <- (usd ^ 2) / umean
# generate a sequence of values for x
x <- seq(0, 1e-8, length = 1000)
# calculate the probability density function (PDF) for the gamma distribution
pdfgamma <- dgamma(x, shape, 1/scale)
# plot
plot(x, pdfgamma, type = "l", lwd = 2, col = "blue",
     xlab = "X", ylab = "PDF", main = "Gamma Distribution")
abline(v = umean, col = "red", lwd = 2, lty = 2)
cal_msci <- msc2time.r(mcmc_msci,
                       u.mean = umean,
                       u.sd = usd,
                       g.mean = gmean, g.sd = gsd)
write.table(cal_msci, "mcmc_comb_cal.txt", quote = FALSE,
            sep = "\t", row.names = FALSE)

# posterior means and 95% HPD
means_msci <- as.data.frame(apply(cal_msci, 2, mean))
postHPD_msci <- as.data.frame(coda::HPDinterval(coda::as.mcmc(cal_msci), ))
summary_msci <- data.frame(parameter = rownames(means_msci),
                           lower95HPD = postHPD_msci$lower,
                           mean = means_msci[, 1],
                           upper95HPD = postHPD_msci$upper)
write.table(summary_msci, "summary_msci.csv", quote = FALSE,
            sep = "\t", row.names = FALSE)

# assemble posterior distribution of trees with divergence times
mcmc <- read.table("mcmc_comb_cal.txt", header = TRUE)
matrix(colnames(mcmc))
mcmc_t <- mcmc[, c(1:11)] # keep only the "t_*" columns
mcmc_t_Ma <- mcmc_t / 1e6 # convert ot Ma

mcmc_t_Ma <- mcmc_t_Ma[sample(nrow(mcmc_t_Ma), 10000), ]

# I use the following lines to add the posterior sample of parameters to the nodes of the tree,
# the code is suitable for this specific species tree, if you want to adjust to your own tree,
# be aware of the order of the parentheses,
# if you have any questions, please send me an email!

tree <- "((som.,att.),((bag.,Sou.),((cha.,Wes.),(Des.,kin.))));"
plot(read.tree(text = tree))

trees_with_times <- matrix(nrow = nrow(mcmc_t_Ma), ncol =  1)
for (i in 1:nrow(trees_with_times)) { # tree topologies
  trees_with_times[i, ] <- tree
}

matrix(colnames(mcmc_t_Ma))
for (i in 1:nrow(trees_with_times)){
  # terminals
  trees_with_times[i] <-  gsub("kin.", paste("kin.:", mcmc_t_Ma[i, 11]), trees_with_times[i])
  trees_with_times[i] <-  gsub("Des.", paste("Des.:", mcmc_t_Ma[i, 11]), trees_with_times[i])
  trees_with_times[i] <-  gsub("Wes.", paste("Wes.:", mcmc_t_Ma[i, 10]), trees_with_times[i])
  trees_with_times[i] <-  gsub("cha.", paste("cha.:", mcmc_t_Ma[i, 10]), trees_with_times[i])
  trees_with_times[i] <-  gsub("Sou.", paste("Sou.:", mcmc_t_Ma[i, 5]), trees_with_times[i])
  trees_with_times[i] <-  gsub("bag.", paste("bag.:", mcmc_t_Ma[i, 5]), trees_with_times[i])
  trees_with_times[i] <-  gsub("som.", paste("som.:", mcmc_t_Ma[i, 2]), trees_with_times[i])
  trees_with_times[i] <-  gsub("att.", paste("att.:", mcmc_t_Ma[i, 2]), trees_with_times[i])
  
  # Node times
  # difference between parent node and mrca node
  # 1st closing parentheses -> (som.,att.)
  trees_with_times[i] <-  str_replace_n(trees_with_times[i], ")", paste("):", mcmc_t_Ma[i, 1] - mcmc_t_Ma[i, 2]), 1)
  # 2nd closing parentheses -> (bag.,Sou.)
  trees_with_times[i] <-  str_replace_n(trees_with_times[i], ")", paste("):", mcmc_t_Ma[i, 4] - mcmc_t_Ma[i, 5]), 2)
  # 3rd closing parentheses -> (cha.,Wes.)
  trees_with_times[i] <-  str_replace_n(trees_with_times[i], ")", paste("):", mcmc_t_Ma[i, 8] - mcmc_t_Ma[i, 10]), 3)
  # 4th closing parentheses -> (Des.,kin.)
  trees_with_times[i] <-  str_replace_n(trees_with_times[i], ")", paste("):", mcmc_t_Ma[i, 8] - mcmc_t_Ma[i, 11]), 4)
  # 5th closing parentheses -> ((cha.,Wes.),(Des.,kin.))
  trees_with_times[i] <-  str_replace_n(trees_with_times[i], ")", paste("):", mcmc_t_Ma[i, 4] - mcmc_t_Ma[i, 8]), 5)
  # 7th closing parentheses -> ((bag.,Sou.),((cha.,Wes.),(Des.,kin.)))
  trees_with_times[i] <-  str_replace_n(trees_with_times[i], ")", paste("):", mcmc_t_Ma[i, 1] - mcmc_t_Ma[i, 4]), 6)
  # 8th closing parenthesis -> ingroup
  trees_with_times[i] <-  str_replace_n(trees_with_times[i], ")", paste("):", mcmc_t_Ma[i, 1] - mcmc_t_Ma[i, 1]), 7)
}

plot(read.tree(text = trees_with_times[[100]]))

# Convert to multiphylo and export
phylo_objects <- lapply(trees_with_times, function(tree) read.tree(text = tree))
multiphylo_object <- as.multiPhylo(phylo_objects)
writeNexus(multiphylo_object, file = "post_dist_caltrees.trees")
# generate maximum clade credibility tree using TreeAnnotator, specifying "median" heights




# TreeMix ----

# NOTE: I opted to not include the TreeMix results in the ms because the results were inconclusive,
# however the commands might come in handy for someone who wants to run the analysis

# library(OptM); library(maditr)
# source("TreeMix_functions.R") # path to required functions for this analysis
# obtained from: "https://github.com/carolindahms/TreeMix"
# setwd("/media/kevin/KEVIN_HDD/software/paquetes_R/BITE/BITE/R/"); sapply(list.files(), source) # source functions from BITE package

# popmap <- read.csv("net_pops.tsv", sep = "\t")
# treemix.from.vcf("populations.snps.maf0.05.vcf", popmap)

#### A. Test migration events ####

# test_linear <- optM("test_migrations/", method = "linear", tsv = NULL) # test m: produces a list of which the $out dataframe suggests optimum m based on multiple linear models
# write.table(test_linear$out, "test_migrations.txt", sep = "\t", quote = FALSE)
# plot_optM(test_linear, method = "linear") # shows changes in log likelihood for different m, and suggests optimum m as 'change points'
# plot 5x5 

# Choose optimum number of m and continue with step 3 in the TreeMix pipeline

#### B. Plot tree ####

# 1. From the final runs, compare tree likelihoods, select tree with highest likelihood, remove duplicates and retain tree(s) with unique topology.
# Adapted from R functions written by Zecca, Labra and Grassi (2019)
# maxLL("final_runs/m2_2m_finalrun_", nt = 10) # first argument is stem of TreeMix output files, nt = number of runs
# shows ML trees and highest likelihood, as well as tree(s) with unique topology. Outputs "TreeLLs.txt" into working directory

# If the number of unique trees is 1, continue with step #2, if n > 1, you might want to create a consensus tree 
# Note that bootstrap and migration values will not be available for consensus tree, thus you could also choose one ML tree 

# cfTrees("final_runs/m2_2m_finalrun_", nt = 10, p = 1, m = "PH85") # m is the method to calculate pairwise distances between unique trees (default is Robinson-Foulds distance)
# p (number from 0.5 to 1) - proportion for a clade to be represented in the consensus tree. Default (1) is strict tree, 0.5 for majority-rule
# plots consensus tree and saves "Consensus.newick" into working directory

# 2. Now plot and save unique tree with highest likelihood
# export PDF lscape 6x6
# treemix.bootstrap("final_runs/m2_2m_finalrun_1", out.file = "tmp", # stem of TreeMix files + run with highest LL and unique topology
#                   phylip.file = "m2_finalconstree.newick",         # consensus tree from the bootstrap procedure generated with PHYLIP
#                   nboot = 1000, fill = TRUE,                       # number of bootstraps
#                   pop.color.file = "data/poporder_col.txt",        # first column is pop name, second the colour name/code (without quotes)
#                   boot.legend.location = "topright")
# pairwise matrix for drift estimates with specified order of pops
# treemix.drift(in.file = "final_runs/m3_3m_finalrun_1", pop.order.color.file = "data/poporder_col.txt") + title("Drift")
# pairwise matrix for residuals
# plot_resid("final_runs/m3_3m_finalrun_1", pop_order = "data/poporder.txt") + title("Residuals")

# 3. If the linear models support multiple optimum admixture models, we can perform a LRT
# change points were detected between 1-5 admx edges
# llik_m1 <- 224.146
# llik_m2 <- 365.365
# llik_m3 <- 433.126
# llik_m4 <- 456.595
# llik_m5 <- 482.076
# 
# pchisq(-2 * (llik_m1 - llik_m2), # -2 * [loglikelihood(nested) - loglikelihood(complex)]; significant p-value = select complex model
#        df = 2, lower.tail = FALSE)                            # 4.670545e-62
# pchisq(-2 * (llik_m2 - llik_m3), df = 2, lower.tail = FALSE)  # 3.730539e-30
# pchisq(-2 * (llik_m3 - llik_m4), df = 2, lower.tail = FALSE)  # 6.420115e-11
# pchisq(-2 * (llik_m4 - llik_m5), df = 2, lower.tail = FALSE)  # 8.58504e-12

#### C. Weights, std. err. and p-values ####

# Output reports mean weight of edge, the jackknife estimate of the weight and standard error (averaged over N independent runs), 
# the least significant p-value recovered  over N runs for each migration event. Set working directory to the final_runs folder.

# GetMigrStats(input_stem = "final_runs/m3_3m_finalrun_", nt = 10) # arguments as above, writes file "MS_and_stats.txt" into current directory

#### Migration support and corrected MS (Optional) ####

# From bootstrap replicates, few other support statistic might be calculated
# The MS is the percentage of times each pair of label sets is present among n bootstrap replicates
# Calculated as: (number of matches/number of bootstrap replicates)*100 from all independent runs in the current working directory
# For the Extended MS (MSE), the number of counts is corrected for multiple matches to avoid over-counting
# Based on R funcions written by Zecca, Labra and Grassi, 2019

# GetPairsOfSets(skipL = 1) # create pairs of sets of populations/taxa from TreeMix output (with treeout.gz extension) in final_runs folder, writes p"PairsOfSets.txt" file
# if you used the flag -noss, set skipL = 2 (default is 1) - the number of lines to skip before reading the tree

# Set working directory to folder with all bootstrap replicates generated with optimum number of m in Step 3
# setwd("~/final_runs/bootstrap")

# Copy PairsOfSets.txt into directory
# GetMigrSupp(skipL = 1) # calculates MS over all bootstrap replicates, writes file "MigrSupp.txt" into current directory
# GetMS_MSe(nmigr = 2, min_n = 2, fixed = "To", skipL = 2) # default input file is "MigrSupp.txt" created with GetMigrSupp(), writes file "MS_MSE.txt" into working directory
# nmigr = number of migrations, fixed = specifies which taxa/species label set of each pair is kept fixed
# fixed = "From" fixes the origin of m; fixed = "To" (default) fixes the destination of the same m 
# min_n = minimum number of taxa/species labels to be included within the unfixed set(s)

# Outputs table with columns 'From' (subset of species below the origin of migration edges), 'To' (the subset of species below the destination of migration edges), Migration Support (MS) and corrected MS with respect to bootstraps (MSE)