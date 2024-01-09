#=
Pkg.add("PhyloNetworks") # add a package
Pkg.add("PhyloPlots") # visualize networks
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("RCall")
Pkg.add("QuartetNetworkGoodnessFit") # goodness of fit test on nets (not used here)
Pkg.add("Gadfly")
Pkg.add("JLD") # save Julia session
=#
using CSV, DataFrames, PhyloNetworks, RCall, PhyloPlots;

cd("phylonetworks")

SNPs2CF = readTableCF("CF_btw_sp.csv")
sppTree = readTopology("(lineomaculatus,(((baguali,south),((west,chacabucoense),(Deseado,kingii))),(som_upt,attenboroughi)));");

net0 = snaq!(sppTree, SNPs2CF, hmax = 0, filename = "btw_sp_net0", runs = 20, seed = 233342);
net0 = readSnaqNetwork("btw_sp_net0.out")
rootatnode!(net0, "lineomaculatus")
plot(net0, :R, usedgelength = true, showedgenumber = true);

net1 = snaq!(net0, SNPs2CF, hmax = 1, filename = "btw_sp_net1", runs = 20, seed = 5634);
net1 = readSnaqNetwork("btw_sp_net1.out")
rootatnode!(net1, "lineomaculatus")
plot(net1, showgamma = true, style = :majortree, showedgenumber = false, arrowlen = 0.2)
# h1 network OK

net2 = snaq!(net1, SNPs2CF, hmax = 2, filename = "btw_sp_net2", runs = 20, seed = 36434);
net2 = readSnaqNetwork("btw_sp_net2.out")
nets2 = readMultiTopology("btw_sp_net2_cp.out")
rootatnode!(nets2[8], "lineomaculatus")
plot(nets2[9], showgamma = true, style = :majortree, showedgenumber = true, arrowlen = 0.2)
# max pseudoL network with 1 reticulation
# the only 2-reticulation network is not compatible with the outgroup rooting

# Given the correct placement of the reticulation in the h1 network inferred in BPP,
# we will optimize branch lengths and inheritance probabilities in this (corrected) network 

net1topo = readTopology("(som_upt,((((kingii,Deseado),((chacabucoense,west),#H10)),(baguali,south)),lineomaculatus),(attenboroughi)#H10);");
rootatnode!(net1topo, "lineomaculatus")
plot(net1topo, style = :majortree, arrowlen = 0.2)
net1par = topologyMaxQPseudolik!(net1topo, SNPs2CF, xtolRel = 1e-10, xtolAbs = 1e-10)
net1par.loglik # pseudo deviance, actually
writeTopology(net1par, "btw_sp_net1opt.out") # export

# pseudo-loglik scores
scores = [net0.loglik, net1par.loglik]
hmax = collect(0:1)
R"plot"(hmax, scores, type = "b", ylab = "network score", xlab = "hmax", col = "blue")