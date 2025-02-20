using Catalyst, Plots, Interpolations, JumpProcesses
using Statistics, Distributions, StatsBase
using CSV, DataFrames, FFTW, DSP
using SymbolicIndexingInterface: parameter_values, state_values
using SciMLStructures: Tunable, replace, replace!
using JLD2, UnPack
using XLSX
using KernelDensity

@parameters kOn, kOff, kteton, ktetoff, kTr, kTrOn, kTrOff, kTrF, kTl, dM, dG;
@variables t;
@species A(t), B(t), DNAoff(t), DNAon(t), DNAonF(t), RNA(t), GFP(t);

# Initial conditions and base parameters
tspan = (0.0, 2000);
u0 = (DNAoff => 1, DNAon => 0, A => 1, B => 1, DNAonF => 0, RNA => 0, GFP => 0); #Changed to approximately level of GFP and RNA at steady state
time = round.(collect(0:215).*0.33334, digits=2);

#p = (kOn => 0.000066*60*60, kOff => 0.000066*60*60,  kTr => 50, kTrOn => 0.00033*60*60, kTrOff => 0.00033*60*60, kTrF => 625, kTl => 475, dM => 1.25, dG => 0.00365)
p = (kOn => 0.00033*60*60, kOff => 0.00066*60*60,  kTr => 15, kTrOn => 0.00025*60*60, kTrOff => 0.00033*60*60, kTrF => 500, kTl => 550, dM => 1.00, dG => 0.00365) 
param_string = join(["$key:$(round(value, digits=2))" for (key, value) in p], ", ")#parameter set for NO FEEDBACK Circuit
ukOn = Normal(0.00033*60*60,0.3*0.00033*60*60) # standard deviation raised to 20%. estimated from Suter et al. paper that used/measured kon and koff values (no reported standard deviation)
ukOff = Normal(0.00066*60*60,0.3*0.00066*60*60)
ukTrOn = Normal(0.00025*60*60,0.3*0.00025*60*60)
ukTrOff = Normal(0.00033*60*60, 0.3*0.00033*60*60)
#udG = Normal(0.00365,0.00126);



#ukTr = Normal(310,150)
# calculate gamma distribution for kTL with mean=550 and std=145
g_scale = 145^2/550  #adjusted from 265^2/500  115^2/475
g_shape = 550/g_scale
ukTl = Gamma(g_shape,g_scale)

kTr_g_scale = 7.5^2/15
kTr_g_shape =15/kTr_g_scale
ukTrg = Gamma(kTr_g_shape,kTr_g_scale)

kTrF_g_scale = 85^2/500
kTrF_g_shape = 500/kTrF_g_scale
ukTrFg = Gamma(kTrF_g_shape,kTrF_g_scale)

kTl_x = 0:01:1000
kTl_y = pdf.(ukTl,kTl_x)
plot(kTl_x, kTl_y, ylabel="PDF",xlabel="kTl")

kTr_x = 0:01:1000
kTr_y = pdf.(ukTrg ,kTr_x)
plot(kTr_x, kTr_y, label="kTrS", ylabel="PDF", xlabel="kTrS (Hr^-1)")

kTrF_x = 0:01:1500
kTrF_y = pdf.(ukTrFg ,kTrF_x)
plot(kTrF_x, kTrF_y, label="kTrF", ylabel="PDF", xlabel="kTrF (Hr^-1)",xticks=0:250:1500)

udG_ln = LogNormal(-5.7,0.42)
dG_x = 0:0.00001:0.25
dG_y = pdf.(udG_ln ,dG_x)
plot(dG_x, dG_y, label="dG", ylabel="PDF", xlabel="dG (Hr^-1)")

function bounded_sample_min0(distribution)
    while true
        sample = rand(distribution)
        if 0 <= sample
            return sample
        end
    end
end
dist_kOn = [];
dist_kOff = [];
dist_kTrOn = [];
dist_kTrOff = [];
dist_kTrg = [];
dist_kTrFg = [];
dist_kTr = [];
dist_dG = [];
for i in 1:1000
    # Generate a new sample
    kOn_new = bounded_sample_min0(ukOn)
    push!(dist_kOn,kOn_new)
    kOff_new = bounded_sample_min0(ukOff)
    push!(dist_kOff,kOff_new)
    kTrOn_new = bounded_sample_min0(ukTrOn)
    push!(dist_kTrOn, kTrOn_new)
    kTrOff_new = bounded_sample_min0(ukTrOff)
    push!(dist_kTrOff, kTrOff_new)
    #kTr_new = rand(ukTr)
    #push!(dist_kTr,kTr_new)  
    kTr_newg = bounded_sample_min0(ukTrg)
    push!(dist_kTrg,kTr_newg)
    kTrF_newg = bounded_sample_min0(ukTrFg)
    push!(dist_kTrFg,kTrF_newg)
    dG_new = bounded_sample_min0(udG_ln)
    push!(dist_dG,dG_new)
end
histogram(dist_kOn,bins=200,label="kOn:$(round(mean(dist_kOn), digits=2))")
histogram(dist_kOff,bins=200,label="kOff:$(round(mean(dist_kOff), digits=2))")
histogram(dist_kTrOn,bins=200,label="kTrOn:$(round(mean(dist_kTrOn), digits=2))")
histogram(dist_kTrOff,bins=200,label="kTrOff:$(round(mean(dist_kTrOff), digits=2))")
#histogram(dist_kTr,bins=200,label="kTr")
histogram(dist_kTrg,bins=200,label="kTrg:$(round(mean(dist_kTrg), digits=2))")
histogram(dist_kTrFg,bins=200,label="kTrgF:$(round(mean(dist_kTrFg), digits=2))")
histogram(dist_dG,bins=200,label="dG:$(round(mean(dist_dG), digits=2))")
# Define maximum value
max_value=[]
min_value=250
max_value = 1000

# Function to sample with a min and maximum value
function bounded_sample(distribution, max_value, min_value)
    while true
        sample = rand(distribution)
        if min_value <= sample <= max_value
            return sample
        end
    end
end
#create ukTl distribution with a maximum value
dist_b_kTl=[];
dist_b_kTr=[];
dist_b_kTrF=[];
for i in 1:1000
    # Generate a new sample
    kTl_new = bounded_sample(ukTl, max_value, min_value)
    kTr_new = bounded_sample(ukTrg,250, 0)
    kTrF_new = bounded_sample(ukTrFg,1250, 250)
    push!(dist_b_kTl,kTl_new)
    push!(dist_b_kTr,kTr_new)
    push!(dist_b_kTrF,kTrF_new)
end
#visuallize the distribution sampled from and the resultant distribution of selected kTl values, representative of what will occur during the simulation loop
#kTl_x = 0:01:1000
#kTl_y = pdf.(ukTl,kTl_x)
#plot(kTl_x, kTl_y)
histogram(dist_b_kTl,bins=200, color=:green,label="kTl:$(round(mean(dist_b_kTl), digits=2))",xlabel="kTl",ylabel="Count")

histogram(dist_b_kTrF,bins=200,color=:blue,label="kTrFast:$(round(mean(dist_b_kTrF), digits=2))",xlabel="kTrFast",ylabel="Count",xlims=(0,1300))
#savefig(joinpath(dir_path, "kTrFast.png"))
histogram(dist_b_kTr,bins=200,color=:red,label="kTrSlow:$(round(mean(dist_b_kTr), digits=2))",xlabel="kTrSlow",ylabel="Count")

histogram(dist_b_kTrF,bins=200,color=:blue,label="kTrFast")
histogram!(dist_b_kTr,bins=200,color=:red,label="kTrSlow")
# Define the reaction network
rxs = [
    (@reaction kOn, DNAoff + A --> DNAon),
    (@reaction kOff, DNAon --> DNAoff + A),
    (@reaction kTr, DNAon --> DNAon + RNA),
    (@reaction kTrOn, DNAoff + B --> DNAonF),
    (@reaction kTrOff, DNAonF --> DNAoff + B),
    (@reaction kTrF, DNAonF --> DNAonF + RNA),
    (@reaction kTl, RNA --> RNA + GFP),
    (@reaction dM, RNA --> 0),
    (@reaction dG, GFP --> 0),
];
# single run before starting loop to check how the "average" trace will look based on parameters    
# Create the ReactionSystem
@named rn = ReactionSystem(rxs, t, [A, B, DNAoff, DNAon, DNAonF, RNA, GFP], [kOn, kOff, kTr, kTrOn, kTrOff, kTrF, kTl, dM, dG]);
rn = complete(rn);

# Define and solve the problem
dprob = DiscreteProblem(rn, u0, tspan, p);
jprob = JumpProblem(rn, dprob, Direct(); save_positions = (false, false));
@time sol = solve(jprob, SSAStepper(); saveat=0.333);

states = unknowns(rn)
params = parameters(rn)
#define steady state section of simulation to where summary statistics will be calculated from
GFP1 = sol[7,:][5250:5465]/1000000
DNAon1 = sol[4,:][5250:5465];
DNAonF1 = sol[5,:][5250:5465];
DNAoff1 = sol[3,:][5250:5465];
#plot GFP and RNA values of the simulation
plot(sol, idxs=7, label="GFP")
plot(GFP1, label=param_string)
plot(sol, idxs=6,label="RNA")
plot(sol, idxs=6,label="RNA",xlims=(1500,1572))
plot(sol, idxs=5,label="DNAonF:$(sum(DNAonF1))",xlims=(1500,1572))
plot(sol, idxs=4,label="DNAon:$(sum(DNAon1))",xlims=(1500,1572))
plot(sol, idxs=3,label="DNAoff:$(sum(DNAoff1))",xlims=(1500,1572))
slope = ((GFP1[end]-GFP1[1])/(time[end]-time[1]))

#setup for Noise Addition
# Estimate the joint KDE
file_path = "C://Users//jsazu//Documents//GitHub//Analysis_of_SGC_Simulations//StreyCats//Broken Feedback Circuit//mean_vs_noise_std.xlsx" #this data is important for adding noise to simulations
data = XLSX.readxlsx(file_path);
sheet = data["Sheet1"];
u = sheet[2:end, 1];
stds = sheet[2:end, 2];

u = convert(Array{Float64}, u);
stds = convert(Array{Float64}, stds);

function sample_cond_dist_std_noise(u, stds, x_fixed; cap_value=77)
    # Cap x_fixed to the maximum allowed value
    capped_x_fixed = min(x_fixed, cap_value)
    data = hcat(u, stds)
    kde_result = kde(data)
    # Create the grid
    ugrid = kde_result.x
    stdgrid = kde_result.y
    # Find the closest index to the capped x value
    u_index = argmin(abs.(ugrid .- capped_x_fixed))
    # Extract the conditional density for the capped x value
    conditional_density = kde_result.density[u_index, :]
    # Normalize the conditional density
    conditional_density /= sum(conditional_density)
    # Sample from the conditional density, ensuring the result is greater than 0
    sampled_y = sample(stdgrid, Weights(conditional_density))
    while sampled_y <= 0
        sampled_y = sample(stdgrid, Weights(conditional_density))
    end

    return sampled_y
end
#Setup for simulation loop
sim_sumstats_list = Vector{Float64}[];
simn_sumstats_list = Vector{Float64}[];
sim_params_list = Vector{Float64}[];
sim_GFP=Vector{Float64}[];
simn_GFP=Vector{Float64}[];
plot();
dir_path = "C:/Users/jsazu/Documents/GitHub/Analysis_of_SGC_Simulations/StreyCats/Broken Feedback Circuit/VaryAct_10_2/" #path to save data and figures

@time for i in 1:1000
    if (i % 100)==0
        println(i)
    end 
    kOn_new = bounded_sample_min0(ukOn)
    kOff_new = bounded_sample_min0(ukOff)
    kTrOn_new = bounded_sample_min0(ukTrOn)
    kTrOff_new = bounded_sample_min0(ukTrOff)
    kTr_new = bounded_sample(ukTrg, 250, 0)
    kTl_new = bounded_sample(ukTl, max_value,min_value)
    kTrF_new = bounded_sample(ukTrFg, 1250, 250)
    dG_new = bounded_sample_min0(udG_ln)
    new_prob = remake(dprob; p = (kOn => kOn_new,
                            kOff => kOff_new,
                            kTr => kTr_new,
                            kTrOn => kTrOn_new,
                            kTrOff => kTrOff_new,
                            kTrF => kTrF_new,
                            kTl => kTl_new,
                            dG => dG_new))
    jprob = JumpProblem(rn, new_prob, Direct(); save_positions = (false, false));
    sol = solve(jprob, SSAStepper(); saveat=0.333) #changed to 0.333 hrs to match measured data
    GFP = sol[7,:][4500:4715]/1000000 #Changed to get values from indexes corresponding to past 1000hrs (steady state) now that there are 3x as many points saved in simulation output (only storing 216 data points) 
    DNAon1 = sol[4,:][4500:4715];
    DNAonF1 = sol[5,:][4500:4715];
    DNAoff1 = sol[3,:][4500:4715];
    # summary statistics of output w/o noise 
    slope = ((GFP[end]-GFP[1])/(time[end]-time[1]))
    k1 = mean(GFP)
    k2 = var(GFP)
    k3 = mean((GFP.-k1).^3)
    k4 = mean((GFP.-k1).^4) .- 3*k2^2
    k5 = mean((GFP.-k1).^5) - 10*k3*k2
    cv = k2/k1
    pw = welch_pgram(GFP)
    ps = pw.power[2:11]
    # Calculate total DNAon, DNAonFast, and DNAoff times
    DNAons = sum(DNAon1)
    DNAonFs = sum(DNAonF1)
    DNAoffs = sum(DNAoff1)
    params = Float64[ kOn_new, kOff_new, kTrOn_new, kTrOff_new, kTr_new, kTrF_new, kTl_new, dG_new]
    plot!(sol, idxs=7, label="GFP $i")
    push!(sim_GFP,GFP)
    push!(sim_params_list,params)
    sim_sumstats = Float64[slope,k1,k2,k3,k4,k5,ps...,cv, DNAons,DNAonFs,DNAoffs]
    push!(sim_sumstats_list,sim_sumstats)

    k1o =()
    k1o = mean(GFP)

    #add noise
    y=()
    y = sample_cond_dist_std_noise(u, stds, k1o)
    noise = y.*randn(216)
    GFPN = GFP+noise
    k1n = mean(GFPN)
    k2n = var(GFPN)
    k3n = mean((GFPN.-k1n).^3)
    k4n = mean((GFPN.-k1n).^4) .- 3*k2n^2
    k5n = mean((GFPN.-k1n).^5) - 10*k3n*k2n
    slopen = ((GFPN[end]-GFPN[1])/(time[end]-time[1]))
    cvn = k2n/k1n
    pwn = welch_pgram(GFPN)
    psn = pwn.power[2:11]
    simn_sumstats = Float64[slopen,k1n,k2n,k3n,k4n,k5n,psn...,cvn]
    push!(simn_sumstats_list,simn_sumstats)
    push!(simn_GFP,GFPN)
end
# plot and save simulation outputs
plot!(legend=:false)
savefig(joinpath(dir_path, "FullSimOutput.png"))
plot(time,sim_GFP,legend=:false,xlabel="Time (hrs)",ylabel="#GFP Molecules (10^6)")
savefig(joinpath(dir_path, "SimOutput.png"))
plot(time,simn_GFP,legend=:false,xlabel="Time (hrs)",ylabel="#GFP Molecules (10^6)")
savefig(joinpath(dir_path, "SimOutputN.png"))
sim_sumstats=[];
simn_sumstats=[];
sim_sumstats = reduce(vcat,transpose.(sim_sumstats_list));
simn_sumstats = reduce(vcat,transpose.(simn_sumstats_list));
params = reduce(vcat,transpose.(sim_params_list));
jldsave(joinpath(dir_path,"noise_2kTr_10_2.jld2"); sim_sumstats,simn_sumstats,params) # save simulation summary stats w/o noise, w noise, and parameters for each trace

sim_GFP_v = reduce(hcat,sim_GFP);
simn_GFP_v = reduce(hcat,simn_GFP)
jldsave(joinpath(dir_path,"noise_2kTr_10_2_sim_GFP.jld2"); sim_GFP_v) # save traces w/o noise
jldsave(joinpath(dir_path,"noise_2kTr_10_2_sinm_GFP.jld2"); simn_GFP_v) # save traces w noise added
#Show Distributions of parameters
histogram(params[:,1], bins=200,xlabel="kOnA",ylabel="counts",label="Mean:$(round(mean(params[:,1]), digits=2))")
savefig(joinpath(dir_path, "kOn.png"))
histogram(params[:,2], bins=200,xlabel="kOffA",ylabel="counts",label="Mean:$(round(mean(params[:,2]), digits=2))")
savefig(joinpath(dir_path, "kOff.png"))
histogram(params[:,3], bins=200,xlabel="kOnB",ylabel="counts",label="Mean:$(round(mean(params[:,3]), digits=2))")
savefig(joinpath(dir_path, "kOnFast.png"))
histogram(params[:,4], bins=200,xlabel="kOffB",ylabel="counts",label="Mean:$(round(mean(params[:,4]), digits=2))")
savefig(joinpath(dir_path, "kOffFast.png"))
histogram(params[:,5], bins=200,xlabel="kTrSlow",ylabel="counts",label="Mean:$(round(mean(params[:,5]), digits=2))")
savefig(joinpath(dir_path, "kTrSlow.png"))
histogram(params[:,6], bins=200,xlabel="kTrF",ylabel="counts",label="Mean:$(round(mean(params[:,6]), digits=2))")
savefig(joinpath(dir_path, "kTrF.png"))
histogram(params[:,7], bins=200,xlabel="kTl",ylabel="counts",label="Mean:$(round(mean(params[:,7]), digits=2))")
savefig(joinpath(dir_path, "kTl.png"))
histogram(params[:,8], bins=200,xlabel="dG",ylabel="counts",label="Mean:$(round(mean(params[:,8]), digits=2))")
savefig(joinpath(dir_path, "dG.png"))

# get summary statistics of experimental data for comparison
@unpack ex_sumstats = jldopen("BrokenFeedbackExperimental.jld2");
# compare summary statistic distributions of experimental and simulated + noise data
plot(ex_sumstats[:,3],ex_sumstats[:,2], seriestype=:scatter, color=:red, xlabel="Variance", ylabel="Mean",label="exp")
plot!(simn_sumstats[:,3],simn_sumstats[:,2], seriestype=:scatter, color=:grey, xlabel="Variance", ylabel="Mean",label="sim")
savefig(joinpath(dir_path, "VarvMean.png"))
plot(ex_sumstats[:,3],ex_sumstats[:,1], seriestype=:scatter, color=:red, xlabel="Variance", ylabel="Slope",label="exp")
plot!(simn_sumstats[:,3],simn_sumstats[:,1], seriestype=:scatter, color=:grey, xlabel="Variance", ylabel="Slope",label="sim")
savefig(joinpath(dir_path, "VarvSlope.png"))
plot(ex_sumstats[:,4],ex_sumstats[:,1], seriestype=:scatter, color=:red, xlabel="Skewness", ylabel="Slope",label="exp")
plot!(simn_sumstats[:,4],simn_sumstats[:,1], seriestype=:scatter, color=:grey, xlabel="Skewness", ylabel="Slope",label="sim")
savefig(joinpath(dir_path, "SkewvSlope.png"))
plot(ex_sumstats[:,5],ex_sumstats[:,1], seriestype=:scatter, color=:red, xlabel="Kurtosis", ylabel="Slope",label="exp")
plot!(simn_sumstats[:,5],simn_sumstats[:,1], seriestype=:scatter, color=:grey, xlabel="Kurtosis", ylabel="Slope",label="sim")
savefig(joinpath(dir_path, "KurtvSlope.png"))
plot(ex_sumstats[:,6],ex_sumstats[:,1], seriestype=:scatter, color=:red, xlabel="k5", ylabel="Slope",label="exp")
plot!(simn_sumstats[:,6],simn_sumstats[:,1], seriestype=:scatter, color=:grey, xlabel="k5", ylabel="Slope",label="sim")
savefig(joinpath(dir_path, "KurtvSlope.png"))
plot(ex_sumstats[:,4],ex_sumstats[:,2], seriestype=:scatter, color=:red, xlabel="Skewness", ylabel="Mean",label="exp")
plot!(simn_sumstats[:,4],simn_sumstats[:,2], seriestype=:scatter, color=:grey, xlabel="Skewness", ylabel="Mean",label="sim")
savefig(joinpath(dir_path, "SkewvMean.png"))
plot(ex_sumstats[:,5],ex_sumstats[:,2], seriestype=:scatter, color=:red, xlabel="Kurtosis", ylabel="Mean",label="exp")
plot!(simn_sumstats[:,5],simn_sumstats[:,2], seriestype=:scatter, color=:grey, xlabel="Kurtosis", ylabel="Mean",label="sim")
savefig(joinpath(dir_path, "KurtvMean.png"))
plot(ex_sumstats[:,6],ex_sumstats[:,2], seriestype=:scatter, color=:red, xlabel="k5", ylabel="Mean",label="exp")
plot!(simn_sumstats[:,6],simn_sumstats[:,2], seriestype=:scatter, color=:grey, xlabel="k5", ylabel="Mean",label="sim")
savefig(joinpath(dir_path, "K5vMean.png"))
plot(ex_sumstats[:,7],ex_sumstats[:,2], seriestype=:scatter, color=:red, xlabel="PS1", ylabel="Mean",label="exp")
plot!(simn_sumstats[:,7],simn_sumstats[:,2], seriestype=:scatter, color=:grey, xlabel="PS1", ylabel="Mean",label="sim")
savefig(joinpath(dir_path, "PS1vMean.png"))
plot(ex_sumstats[:,7],ex_sumstats[:,3], seriestype=:scatter, color=:red, xlabel="PS1", ylabel="Variance",label="exp")
plot!(simn_sumstats[:,7],simn_sumstats[:,3], seriestype=:scatter, color=:grey, xlabel="PS1", ylabel="Variance",label="sim")
savefig(joinpath(dir_path, "PS1vVar.png"))
plot(params[:,8],simn_sumstats[:,2], seriestype=:scatter, xlabel="degG", ylabel="Mean")
savefig(joinpath(dir_path, "degGvMean.png"))
plot(params[:,7],simn_sumstats[:,2], seriestype=:scatter, xlabel="kTl", ylabel="Mean")
savefig(joinpath(dir_path, "kTlvMean.png"))
plot(params[:,6],simn_sumstats[:,2], seriestype=:scatter, xlabel="kTrF", ylabel="Mean")
savefig(joinpath(dir_path, "kTrFvMean.png"))
plot(params[:,5],simn_sumstats[:,2], seriestype=:scatter, xlabel="kTr", ylabel="Mean")
savefig(joinpath(dir_path, "kTrvMean.png"))
plot(params[:,4],simn_sumstats[:,2], seriestype=:scatter, xlabel="kFOff", ylabel="Mean")
savefig(joinpath(dir_path, "kFOffvMean.png"))
plot(params[:,3],simn_sumstats[:,2], seriestype=:scatter, xlabel="kFon", ylabel="Mean")
savefig(joinpath(dir_path, "kFonvMean.png"))
plot(params[:,2],simn_sumstats[:,2], seriestype=:scatter, xlabel="kOn", ylabel="Mean")
savefig(joinpath(dir_path, "konvMean.png"))
plot(params[:,1],simn_sumstats[:,2], seriestype=:scatter, xlabel="kOff", ylabel="Mean")
savefig(joinpath(dir_path, "koffvMean.png"))
plot(sim_sumstats[:,18],simn_sumstats[:,2], seriestype=:scatter, xlabel="Slow On time", ylabel="Mean")
savefig(joinpath(dir_path, "slowontimevMean.png"))
plot(sim_sumstats[:,19],simn_sumstats[:,2], seriestype=:scatter, xlabel="Fast On Time", ylabel="Mean")
savefig(joinpath(dir_path, "fastontimevMean.png"))
plot(sim_sumstats[:,20],simn_sumstats[:,2], seriestype=:scatter, xlabel="Off Time", ylabel="Mean")
savefig(joinpath(dir_path, "offtimevMean.png"))
plot(params[:,5],params[:,6], seriestype=:scatter, xlabel="kTr", ylabel="kTrF")
savefig(joinpath(dir_path, "kTrvkTrF.png"))
plot(params[:,6],params[:,7], seriestype=:scatter, xlabel="kTrF", ylabel="kTl")
savefig(joinpath(dir_path, "kTrFvkTl.png"))
plot(params[:,5],params[:,7], seriestype=:scatter, xlabel="kTr", ylabel="kTl")
savefig(joinpath(dir_path, "kTrvkTl.png"))

# Visually compare the summary statistics of the simulated + noise time traces and the experimental time traces
bins=200
# slope
min_edge = min(minimum(simn_sumstats[:,1]), minimum(ex_sumstats[:,1]))
max_edge = max(maximum(simn_sumstats[:,1]), maximum(ex_sumstats[:,1]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,1], normalize=:probability, color=:red, bins=bin_edges, alpha=0.5,label="Exp Slope")
#histogram!(sim_sumstats[:,1],bins=6000,color=:grey,alpha=0.5,label="Sim Temp Means")
histogram!(simn_sumstats[:,1], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Slope", ylabel="Probability", xlabel="Slope of GFP # (10^6)")
savefig(joinpath(dir_path, "SlopeHistExp_vs_Sim_Broken.png"))
#Mean
min_edge = min(minimum(simn_sumstats[:,2]), minimum(ex_sumstats[:,2]))
max_edge = max(maximum(simn_sumstats[:,2]), maximum(ex_sumstats[:,2]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,2], color=:red, normalize=:probability,  bins=bin_edges, alpha=0.5, label="Exp Temp Mean")
histogram!(simn_sumstats[:,2], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Temp Mean", ylabel="Probability", xlabel="Temporal Mean of GFP # (10^6)") #,ylim=(0,0.125),xlim=(0,20)
savefig(joinpath(dir_path, "MeanHistExp_vs_Sim_Broken.png"))
# Variance
min_edge = min(minimum(simn_sumstats[:,3]), minimum(ex_sumstats[:,3]))
max_edge = max(maximum(simn_sumstats[:,3]), maximum(ex_sumstats[:,3]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,3], color=:red, normalize=:probability,  bins=bin_edges, alpha=0.5, label="Exp Temp Vars")
histogram!(simn_sumstats[:,3], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Temp Vars", ylabel="Probability", xlabel="Temporal Variance of GFP # (10^6)") #,ylim=(0,0.125),xlim=(0,20)
savefig(joinpath(dir_path, "VarHistExp_vs_Sim_Broken.png.png"))
min_edge = min(minimum(simn_sumstats[:,3]), minimum(ex_sumstats[:,3]))
max_edge = max(maximum(simn_sumstats[:,3]), maximum(ex_sumstats[:,3]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,3], color=:red, normalize=:probability,  bins=bin_edges, alpha=0.5, label="Exp Temp Vars")
histogram!(simn_sumstats[:,3], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Temp Vars", ylabel="Probability", xlabel="Temporal Variance of GFP # (10^6)",ylim=(0,0.125),xlim=(0,20))
savefig(joinpath(dir_path, "VarzHistExp_vs_Sim_Broken.png"))
# CV
min_edge = min(minimum(simn_sumstats[:,17]), minimum(ex_sumstats[:,17]))
max_edge = max(maximum(simn_sumstats[:,17]), maximum(ex_sumstats[:,17]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,17], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp Temp CV")
histogram!(simn_sumstats[:,17], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Temp Vars", ylabel="Probability", xlabel="Temporal CV of GFP # (10^6)")
savefig(joinpath(dir_path, "CVExp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,17]), minimum(ex_sumstats[:,17]))
max_edge = max(maximum(simn_sumstats[:,17]), maximum(ex_sumstats[:,17]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,17], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp Temp CV")
histogram!(simn_sumstats[:,17], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Temp Vars", ylabel="Probability", xlabel="Temporal CV of GFP # (10^6)",ylim=(0,0.175),xlim=(0,0.5))
savefig(joinpath(dir_path, "CVzExp_vs_Sim_Broken.png"))
#k3 - skewness
min_edge = min(minimum(simn_sumstats[:,4]), minimum(ex_sumstats[:,4]))
max_edge = max(maximum(simn_sumstats[:,4]), maximum(ex_sumstats[:,4]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,4], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k3")
histogram!(simn_sumstats[:,4], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k3", ylabel="Probability", xlabel="k3") #, ylim=(0,0.20), xlim=(-30,50)
savefig(joinpath(dir_path, "SkewExp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,4]), minimum(ex_sumstats[:,4]))
max_edge = max(maximum(simn_sumstats[:,4]), maximum(ex_sumstats[:,4]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,4], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k3")
histogram!(simn_sumstats[:,4], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k3", ylabel="Probability", xlabel="k3", ylim=(0,0.20), xlim=(-30,50))
savefig(joinpath(dir_path, "SkewzExp_vs_Sim_Broken.png"))
#k4 - Kurtosis
min_edge = min(minimum(simn_sumstats[:,5]), minimum(ex_sumstats[:,5]))
max_edge = max(maximum(simn_sumstats[:,5]), maximum(ex_sumstats[:,5]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,5],  color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k4")
histogram!(simn_sumstats[:,5], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k4", ylabel="Probability", xlabel="k4") #, ylim=(0,0.20), xlim=(-500,250)
savefig(joinpath(dir_path, "KurtExp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,5]), minimum(ex_sumstats[:,5]))
max_edge = max(maximum(simn_sumstats[:,5]), maximum(ex_sumstats[:,5]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,5],  color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k4")
histogram!(simn_sumstats[:,5], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k4", ylabel="Probability", xlabel="k4", ylim=(0,0.20), xlim=(-500,250))
savefig(joinpath(dir_path, "KurtzExp_vs_Sim_Broken.png"))
#k5
min_edge = min(minimum(simn_sumstats[:,6]), minimum(ex_sumstats[:,6]))
max_edge = max(maximum(simn_sumstats[:,6]), maximum(ex_sumstats[:,6]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,6],  color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k5")
histogram!(simn_sumstats[:,6], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k5", ylabel="Probability", xlabel="k5") #, ylim=(0,0.20),xlim=(-10000,5000)
savefig(joinpath(dir_path, "K5Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,6]), minimum(ex_sumstats[:,6]))
max_edge = max(maximum(simn_sumstats[:,6]), maximum(ex_sumstats[:,6]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,6],  color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k5")
histogram!(simn_sumstats[:,6], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k5", ylabel="Probability", xlabel="k5", ylim=(0,0.20),xlim=(-10000,5000))
savefig(joinpath(dir_path, "K5zExp_vs_Sim_Broken.png"))
#power spectrum
min_edge = min(minimum(simn_sumstats[:,7]), minimum(ex_sumstats[:,7]))
max_edge = max(maximum(simn_sumstats[:,7]), maximum(ex_sumstats[:,7]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,7], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k6")
histogram!(simn_sumstats[:,7], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k6", ylabel="Probability", xlabel="k6") #, ylim=(0,0.085),xlim=(0,10)
savefig(joinpath(dir_path, "PS1Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,7]), minimum(ex_sumstats[:,7]))
max_edge = max(maximum(simn_sumstats[:,7]), maximum(ex_sumstats[:,7]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,7], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k6")
histogram!(simn_sumstats[:,7], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k6", ylabel="Probability", xlabel="k6", ylim=(0,0.09),xlim=(0,15))
savefig(joinpath(dir_path, "PS1zExp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,8]), minimum(ex_sumstats[:,8]))
max_edge = max(maximum(simn_sumstats[:,8]), maximum(ex_sumstats[:,8]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,8], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k7")
histogram!(simn_sumstats[:,8], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k7", ylabel="Probability", xlabel="k7") #, ylim=(0,0.1),xlim=(0,3)
savefig(joinpath(dir_path, "PS2Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,8]), minimum(ex_sumstats[:,8]))
max_edge = max(maximum(simn_sumstats[:,8]), maximum(ex_sumstats[:,8]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,8], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k7")
histogram!(simn_sumstats[:,8], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k7", ylabel="Probability", xlabel="k7", ylim=(0,0.11),xlim=(0,6))
savefig(joinpath(dir_path, "PS2zExp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,9]), minimum(ex_sumstats[:,9]))
max_edge = max(maximum(simn_sumstats[:,9]), maximum(ex_sumstats[:,9]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,9], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k8")
histogram!(simn_sumstats[:,9], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k8", ylabel="Probability", xlabel="k8") #, ylim=(0,0.1),xlim=(0,2.25)
savefig(joinpath(dir_path, "PS3Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,9]), minimum(ex_sumstats[:,9]))
max_edge = max(maximum(simn_sumstats[:,9]), maximum(ex_sumstats[:,9]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,9], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k8")
histogram!(simn_sumstats[:,9], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k8", ylabel="Probability", xlabel="k8", ylim=(0,0.5),xlim=(0,8))
savefig(joinpath(dir_path, "PS3zExp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,10]), minimum(ex_sumstats[:,10]))
max_edge = max(maximum(simn_sumstats[:,10]), maximum(ex_sumstats[:,10]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,10], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k9")
histogram!(simn_sumstats[:,10], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k9", ylabel="Probability", xlabel="k9") #, ylim=(0,0.1),xlim=(0,1.75)
savefig(joinpath(dir_path, "PS4Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,10]), minimum(ex_sumstats[:,10]))
max_edge = max(maximum(simn_sumstats[:,10]), maximum(ex_sumstats[:,10]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,10], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k9")
histogram!(simn_sumstats[:,10], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k9", ylabel="Probability", xlabel="k9", ylim=(0,0.8),xlim=(0,8))
savefig(joinpath(dir_path, "PS4zExp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,11]), minimum(ex_sumstats[:,11]))
max_edge = max(maximum(simn_sumstats[:,11]), maximum(ex_sumstats[:,11]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,11], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k10")
histogram!(simn_sumstats[:,11], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k10", ylabel="Probability", xlabel="k10") #, ylim=(0,0.1),xlim=(0,1.75)
savefig(joinpath(dir_path, "PS5Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,11]), minimum(ex_sumstats[:,11]))
max_edge = max(maximum(simn_sumstats[:,11]), maximum(ex_sumstats[:,11]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,11], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k10")
histogram!(simn_sumstats[:,11], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k10", ylabel="Probability", xlabel="k10", ylim=(0,0.85),xlim=(0,8))
savefig(joinpath(dir_path, "PS5zExp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,12]), minimum(ex_sumstats[:,12]))
max_edge = max(maximum(simn_sumstats[:,12]), maximum(ex_sumstats[:,12]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,12], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k11")
histogram!(simn_sumstats[:,12], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k11", ylabel="Probability", xlabel="k11") #, ylim=(0,0.225),xlim=(0,3)
savefig(joinpath(dir_path, "PS6Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,12]), minimum(ex_sumstats[:,12]))
max_edge = max(maximum(simn_sumstats[:,12]), maximum(ex_sumstats[:,12]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,12], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k11")
histogram!(simn_sumstats[:,12], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k11", ylabel="Probability", xlabel="k11", ylim=(0,0.8),xlim=(0,8))
savefig(joinpath(dir_path, "PS6zExp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,13]), minimum(ex_sumstats[:,13]))
max_edge = max(maximum(simn_sumstats[:,13]), maximum(ex_sumstats[:,13]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,13], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k12")
histogram!(simn_sumstats[:,13], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k12", ylabel="Probability", xlabel="k12") #, ylim=(0,0.175),xlim=(0,1.25)
savefig(joinpath(dir_path, "PS7Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,13]), minimum(ex_sumstats[:,13]))
max_edge = max(maximum(simn_sumstats[:,13]), maximum(ex_sumstats[:,13]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,13], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k12")
histogram!(simn_sumstats[:,13], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k12", ylabel="Probability", xlabel="k12", ylim=(0,0.9),xlim=(0,8))
savefig(joinpath(dir_path, "PS7zExp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,14]), minimum(ex_sumstats[:,14]))
max_edge = max(maximum(simn_sumstats[:,14]), maximum(ex_sumstats[:,14]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,14], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k13")
histogram!(simn_sumstats[:,14], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k13", ylabel="Probability", xlabel="k13") #, ylim=(0,0.175),xlim=(0,1.25)
savefig(joinpath(dir_path, "PS8Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,14]), minimum(ex_sumstats[:,14]))
max_edge = max(maximum(simn_sumstats[:,14]), maximum(ex_sumstats[:,14]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,14], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k13")
histogram!(simn_sumstats[:,14], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k13", ylabel="Probability", xlabel="k13", ylim=(0,0.8),xlim=(0,8))
savefig(joinpath(dir_path, "PS8zExp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,15]), minimum(ex_sumstats[:,15]))
max_edge = max(maximum(simn_sumstats[:,15]), maximum(ex_sumstats[:,15]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,15], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k14")
histogram!(simn_sumstats[:,15], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k14", ylabel="Probability", xlabel="k14") #, ylim=(0,0.185),xlim=(0,1.85)
savefig(joinpath(dir_path, "PS9Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,15]), minimum(ex_sumstats[:,15]))
max_edge = max(maximum(simn_sumstats[:,15]), maximum(ex_sumstats[:,15]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,15], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k14")
histogram!(simn_sumstats[:,15], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k14", ylabel="Probability", xlabel="k14", ylim=(0,0.8),xlim=(0,8))
savefig(joinpath(dir_path, "PS9zExp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,16]), minimum(ex_sumstats[:,16]))
max_edge = max(maximum(simn_sumstats[:,16]), maximum(ex_sumstats[:,16]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,16], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label=" Exp k15")
histogram!(simn_sumstats[:,16], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k15", ylabel="Probability", xlabel="k15") #, ylim=(0,0.225),xlim=(0,2)
savefig(joinpath(dir_path, "PS10Exp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,16]), minimum(ex_sumstats[:,16]))
max_edge = max(maximum(simn_sumstats[:,16]), maximum(ex_sumstats[:,16]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,16], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label=" Exp k15")
histogram!(simn_sumstats[:,16], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k15", ylabel="Probability", xlabel="k15", ylim=(0,0.8),xlim=(0,8))
savefig(joinpath(dir_path, "PS10zExp_vs_Sim_Broken.png"))