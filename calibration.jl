using CSV
using DataFrames
using DelimitedFiles
using PhyloNetworks
using Pkg
using RCall
using Statistics
using StatsBase

cd("phylonetworks/")

net0 = rootatnode!(readSnaqNetwork("btw_sp_net0.out"), "lineomaculatus")
net1 = rootatnode!(readTopology("btw_sp_net1opt.out"), "lineomaculatus")
# check the rooting
plot(net0, :R)
plot(net1)

taxa = tipLabels(net1)
ntips = length(taxa) # 9
outgroup = Set(["lineomaculatus"])
ingroup  = setdiff(tipLabels(net1), outgroup)

#########################
## Network calibration ##
#########################

# An option to calibrate a network (time-consistent) is to estimate pariwise distances between 
# terminals from estimated gene trees. Branch lengths of such trees must represent expected substitutions.
# However, short DNA fragments (such as RAD loci) have a low information content to estimate branch lengths.
# In these cases, it might be best to estimate pairwise genetic distances directly from the SNPs matrix:

# Here, I estimated distances from a VCF file
# the procedure is detailed in the R script `phylogeography.R` 

# type `$` in the Julia REPL to access R
R"library(vcfR); library(adegenet)
vcf_inds <- read.vcfR('../data/genomics/populations.snps.maf0.01.net.vcf')
gen_inds <- vcfR2genind(vcf_inds)
ploidy(gen_inds) <- 2
gen_pops <- genind2genpop(gen_inds, pop = factor(read.table('../data/sheets/net_pops.tsv')[, 2]))
avdist <-  as.matrix(dist.genpop(gen_pops))
# match the row and column names with the tip order of the networks
# use `$` to retrieve a Julia variable into R
avdist <- avdist[match($taxa, rownames(avdist)), match($taxa, colnames(avdist))]"

# back in Julia, do `rcopy` to retrieve an R variable into Julia
avdist = Matrix(rcopy(R"avdist"))
avdist |> maximum # max distance: ≈ 0.134 substitutions/site, averaged across genes


# net1

net1_cp = deepcopy(net1)
# modify branch lengths for good starting values: max distance / 9 taxa
for e in net1_cp.edge
    e.length = (e.hybrid ? 0.0 : (avdist |> maximum) / length(taxa))
end
using Random; Random.seed!(8732) # for reproducibility: you can pick any other seed
fmin, xmin, ret = calibrateFromPairwiseDistances!(net1_cp, avdist, taxa, forceMinorLength0 = false, NLoptMethod = :LD_MMA);
fmin, xmin, ret = calibrateFromPairwiseDistances!(net1_cp, avdist, taxa, forceMinorLength0 = false, NLoptMethod = :LN_COBYLA);
fmin, xmin, ret = calibrateFromPairwiseDistances!(net1_cp, avdist, taxa, forceMinorLength0 = false, NLoptMethod = :LD_MMA);
sort!([e.length for e in net1_cp.edge])

# replace negative values with a very small value
for e in net1_cp.edge
    if e.length < 0 e.length = 0.01; end
end
plot(net1_cp, useedgelength = true, showedgenumber = true);

# normalize the network height
rootage = getNodeAges(net1_cp)[1] # ~ 0.002, ~ 0.136 con gtreePairDist_vcf
for e in net1_cp.edge
    e.length /= rootage
end
getNodeAges(net1_cp) |> maximum # 1.0, sanity check
plot(net1_cp, useedgelength = true, style = :majortree, arrowlen = 0.2);


# You can do the same for the tree (net0) if you want

# save the calibrated networks
writeTopology(net1_cp, "net1_cal_distsnps.phy")

# Plot the networks

R"pdf"("plots_net1.pdf", width = 7, height = 5);
R"par"(mar = [0, 0, 0, 0]);
PhyloPlots.plot(net1, useedgelength = true, style = :majortree, arrowlen = 0.2, showedgelength = true, showgamma = true);
PhyloPlots.plot(net1_cp, useedgelength = true, style = :majortree, arrowlen = 0.2, showedgelength = true, showgamma = true);
R"dev.off()";
