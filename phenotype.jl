using CSV
using DataFrames
using PhyloNetworks
using PhyloPlots
using RCall
using StatsBase
using StatsModels

cd("phylonetworks/")

net1 = readTopology("net1_cal_distsnps.phy");
tree  = majorTree(net1)
plot(tree, useedgelength = true, showedgelength = true)

traitsM = CSV.read("data/phenotype/phenoevol_m_allom.csv", DataFrame)
traitsM = traitsM[:, [:sp, :SVL, :AL, :TFL, :HH, :HW, :RND, :alt, :lat, :long]]
traitsF = CSV.read("data/phenotype/phenoevol_f_allom.csv", DataFrame)
traitsF = traitsF[:, [:sp, :SVL, :AL, :TFL, :HH, :HW, :RND, :alt, :lat, :long]]
traits = CSV.read("data/phenotype/phenoevol_allom.csv", DataFrame)
traits = traits[:, [:sp, :SVL, :AL, :TFL, :HH, :HW, :RND, :alt, :lat, :long]]
geom = CSV.read("data/phenotype/phenoevol_gm_pcscores.csv", DataFrame)

subset!(traitsM, :sp => ByRow(in(tipLabels(net1))));
subset!(traitsF, :sp => ByRow(in(tipLabels(net1))));
subset!(traits, :sp => ByRow(in(tipLabels(net1))));
subset!(geom, :sp => ByRow(in(tipLabels(net1))));

# data processing 
traits_spp = combine(groupby(traitsF, "sp"),
    "SVL"	=> mean	=> "svl",
	"SVL"	=> std	=> "svl_sd", nrow => :svl_n,
	"AL"	=> mean	=> "al",
	"AL"	=> std	=> "al_sd", nrow => :al_n,
	"TFL"	=> mean	=> "tfl",
	"TFL"	=> std	=> "tfl_sd", nrow => :tfl_n,
    "HH"	=> mean	=> "hh",
	"HH"	=> std	=> "hh_sd", nrow => :hh_n,
    "HW"	=> mean	=> "hw",
	"HW"	=> std	=> "hw_sd", nrow => :hw_n,
	"RND"	=> mean	=> "rnd",
	"RND"	=> std	=> "rnd_sd", nrow => :rnd_n,
	"alt"	=> mean	=> "alt",
    "lat"	=> mean	=> "lat",
	"long"	=> mean	=> "long")
	# repeat with the male and female subsets
geom_spp = combine(groupby(geom, "sp"),
	"Comp1" => mean => "pc1",
	"Comp1" => std => "pc1_sd", nrow => :pc1_n,
	"Comp2" => mean => "pc2",
	"Comp2" => std => "pc2_sd", nrow => :pc2_n,
	"Comp3" => mean => "pc3",
	"Comp3" => std => "pc3_sd", nrow => :pc3_n,
	"Comp4" => mean => "pc4",
	"Comp4" => std => "pc4_sd", nrow => :pc4_n,
	"Comp5" => mean => "pc5",
	"Comp5" => std => "pc5_sd", nrow => :pc5_n,
	"Comp6" => mean => "pc6",
	"Comp6" => std => "pc6_sd", nrow => :pc6_n,
	"Comp7" => mean => "pc7",
	"Comp7" => std => "pc7_sd", nrow => :pc7_n,
	"Comp8" => mean => "pc8",
	"Comp8" => std => "pc8_sd", nrow => :pc8_n,
	"Comp9" => mean => "pc9",
	"Comp9" => std => "pc9_sd", nrow => :pc9_n,
	"Comp10" => mean => "pc10",
	"Comp10" => std => "pc10_sd", nrow => :pc10_n,
	"Comp11" => mean => "pc11",
	"Comp11" => std => "pc11_sd", nrow => :pc11_n)

CSV.write("pheno_allom_proc.csv", traits_spp) # repeat with the male and female subsets
CSV.write("pheno_geom_proc.csv", geom_spp)

traitsM = CSV.read("pheno_f_allom_proc.csv", DataFrame)
traitsF = CSV.read("pheno_m_allom_proc.csv", DataFrame)
traits = CSV.read("pheno_allom_proc.csv", DataFrame)
geom = CSV.read("pheno_geom_proc.csv", DataFrame)

############################
## Network vs tree models ##
############################

# without predictors
fit_tree = phylolm(@formula(pc11 ~ 1), geom, tree, tipnames = :sp, withinspecies_var = true, y_mean_std = true, reml = false)
fit_net1 = phylolm(@formula(pc11 ~ 1), geom, net1, tipnames = :sp, withinspecies_var = true, y_mean_std = true, reml = false)
DataFrame(logLik = [loglikelihood(fit_tree), loglikelihood(fit_net1)],
		  AICc = [aicc(fit_tree), aicc(fit_net1)],
		  deltaAICc = [aicc(fit_tree) - aicc(fit_tree),
					  aicc(fit_tree) - aicc(fit_net1)])

#############################
## Transgressive evolution ##
#############################

## h1 network
df_shift = regressorHybrid(net1)
# this table shows the branches involved in introgression (shift_* columns),
# column tipNames shows the terminals affected by the introgression

# here you can include correlated spatial effects (alt, lat, long), or not
df = innerjoin(select(geom, :sp, r"pc11",
#					  :alt,
#					  :lat,
#					  :long
					),
			select(df_shift, Not(:sum)), # excludes the 'sum' column
			on = :sp => :tipNames) # join our data with shift predictors

fit_noshift = phylolm(@formula(pc11 ~ 1), df, net1, tipnames = :sp, withinspecies_var = true, y_mean_std = true, reml = false)
# using `df` for both, so the no-shift model is recognized as nested in the shifted model
fit_sh3 = phylolm(@formula(pc11 ~ shift_3), df, net1, tipnames = :sp, withinspecies_var = true, y_mean_std = true, reml = false) # attenboroughi
coeftable(fit_sh3) # 95% confidence interval for shift includes 0?
DataFrame(logLik = [loglikelihood(fit_noshift), loglikelihood(fit_sh3)],
			AICc = [aicc(fit_noshift), aicc(fit_sh3)],
			deltaAICc = [aicc(fit_noshift) - aicc(fit_noshift),
					aicc(fit_noshift) - aicc(fit_sh3)])
# we could run a likelihood ratio test (lrtest) to compare the two models,
# this is only possible with ML optimization
lrtest(fit_noshift, fit_sh3)