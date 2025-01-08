using Catalyst, Plots, Interpolations, JumpProcesses
using Statistics, Distributions, StatsBase
using CSV, DataFrames, FFTW, DSP
using SymbolicIndexingInterface: parameter_values, state_values
using SciMLStructures: Tunable, replace, replace!
using JLD2, UnPack
using XLSX

@parameters kOn, kOff, kteton, ktetoff, kTr, kTrOn, kTrOff, kTrF, kTl, dM, dG;
@variables t;
@species A(t), B(t), DNAoff(t), DNAon(t), DNAonF(t), RNA(t), GFP(t);

# Initial conditions and base parameters
tspan = (0.0, 2000);
u0 = (DNAoff => 1, DNAon => 0, A => 1, B => 1, DNAonF => 0, RNA => 0, GFP => 0); #Changed to approximately level of GFP and RNA at steady state
time = round.(collect(0:215).*0.33334, digits=2);

p = (kOn => 0.0033*60*60, kOff => 0.00033*60*60,  kTr => 50, kTrOn => 0.0033*60*60, kTrOff => 0.00033*60*60, kTrF => 500, kTl => 425, dM => 1.25, dG => 0.00365) 
param_string = join(["$key:$(round(value, digits=2))" for (key, value) in p], ", ")#parameter set for NO FEEDBACK Circuit
ukOn = Normal(0.0033*60*60,0.2*0.0033*60*60) # standard deviation raised to 20%. estimated from Suter et al. paper that used/measured kon and koff values (no reported standard deviation)
ukOff = Normal(0.00033*60*60,0.2*0.00033*60*60)
ukTrOn = Normal(0.0033*60*60,0.2*0.0033*60*60)
ukTrOff = Normal(0.00033*60*60, 0.2*0.00033*60*60)
#ukTr = Normal(310,150)
# calculate gamma distribution for kTL with mean=310 and std=265
g_scale = 80^2/425  #adjusted from 265^2/500
g_shape = 425/g_scale
ukTl = Gamma(g_shape,g_scale)

kTr_g_scale = 47.5^2/70
kTr_g_shape =70/kTr_g_scale
ukTrg = Gamma(kTr_g_shape,kTr_g_scale)

kTrF_g_scale = 80^2/480
kTrF_g_shape = 480/kTrF_g_scale
ukTrFg = Gamma(kTrF_g_shape,kTrF_g_scale)

kTr_x = 0:01:1000
kTr_y = pdf.(ukTrg ,kTr_x)
plot(kTr_x, kTr_y)

kTrF_x = 0:01:1000
kTrF_y = pdf.(ukTrFg ,kTrF_x)
plot(kTrF_x, kTrF_y)

dist_kOn = [];
dist_kOff = [];
dist_kTrg = [];
dist_kTrFg = [];
dist_kTr = [];
for i in 1:1000
    # Generate a new sample
    kOn_new = rand(ukOn)
    push!(dist_kOn,kOn_new)
    kOff_new = rand(ukOff)
    push!(dist_kOff,kOff_new) 
    kTr_new = rand(ukTr)
    push!(dist_kTr,kTr_new)  
    kTr_newg = rand(ukTrg)
    push!(dist_kTrg,kTr_newg)
    kTrF_newg = rand(ukTrFg)
    push!(dist_kTrFg,kTrF_newg)
end
histogram(dist_kOn,bins=200,label="kOn")
histogram(dist_kOff,bins=200,label="kOff")
histogram(dist_kTr,bins=200,label="kTr")
histogram(dist_kTrg,bins=200,label="kTrg")
histogram(dist_kTrFg,bins=200,label="kTrgF")
# Define maximum value
max_value=[]
min_value=100
max_value = 750

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
    kTr_new = bounded_sample(ukTrg,100, 0)
    kTrF_new = bounded_sample(ukTrFg,750, 200)
    push!(dist_b_kTl,kTl_new)
    push!(dist_b_kTr,kTr_new)
    push!(dist_b_kTrF,kTrF_new)
end
#visuallize the distribution sampled from and the resultant distribution of selected kTl values, representative of what will occur during the simulation loop
kTl_x = 0:01:1000
kTl_y = pdf.(ukTl,kTl_x)
plot(kTl_x, kTl_y)
histogram(dist_b_kTl,bins=200)

histogram(dist_b_kTrF,bins=200,color=:blue,label="kTrFast")
histogram(dist_b_kTr,bins=200,color=:red,label="kTrSlow")

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
GFP1 = sol[7,:][4500:4715]/1000000
#plot GFP and RNA values of the simulation
plot(sol, idxs=7, label="GFP")
plot(GFP1, label=param_string)
plot(sol, idxs=6,label="RNA")
plot(sol, idxs=5,label="DNAonF",xlims=(1500,1572))
plot(sol, idxs=4,label="DNAon",xlims=(1500,1572))
plot(sol, idxs=3,label="DNAoff",xlims=(1500,1572))
slope = ((GFP1[end]-GFP1[1])/(time[end]-time[1]))

#Setup for simulation loop
sim_sumstats_list = Vector{Float64}[]
simn_sumstats_list = Vector{Float64}[]
sim_params_list = Vector{Float64}[]
sim_GFP=Vector{Float64}[]
simn_GFP=Vector{Float64}[]
plot();
@time for i in 1:1000
    if (i % 100)==0
        println(i)
    end 
    kOn_new = rand(ukOn)
    kOff_new = rand(ukOff)
    kTrOn_new = rand(ukTrOn)
    kTrOff_new = rand(ukTrOff)
    kTr_new = bounded_sample(ukTrg, 100, 15)
    kTl_new = bounded_sample(ukTl, max_value,min_value)
    kTrF_new = bounded_sample(ukTrFg, 750, 200)
    new_prob = remake(dprob; p = (kOn => kOn_new,
                            kOff => kOff_new,
                            kTr => kTr_new,
                            kTrF => kTrF_new,
                            kTl => kTl_new))
    jprob = JumpProblem(rn, new_prob, Direct(); save_positions = (false, false));
    sol = solve(jprob, SSAStepper(); saveat=0.333) #changed to 0.333 hrs to match measured data
    GFP = sol[7,:][4500:4715]/1000000 #Changed to get values from indexes corresponding to past 1000hrs (steady state) now that there are 3x as many points saved in simulation output (only storing 216 data points) 
    
    slope = ((GFP[end]-GFP[1])/(time[end]-time[1]))
    k1 = mean(GFP)
    k2 = var(GFP)
    k3 = mean((GFP.-k1).^3)
    k4 = mean((GFP.-k1).^4) .- 3*k2^2
    k5 = mean((GFP.-k1).^5) - 10*k3*k2
    cv = k2/k1
    pw = welch_pgram(GFP)
    ps = pw.power[2:11]
    params = Float64[ kOn_new, kOff_new, kTr_new, kTl_new]
    plot!(sol, idxs=7, label="GFP $i")
    push!(sim_GFP,GFP)
    push!(sim_params_list,params)
    sim_sumstats = Float64[slope,k1,k2,k3,k4,k5,ps...,cv]
    push!(sim_sumstats_list,sim_sumstats)

    k1o =()
    k1o = mean(GFP)

    #add noise
    y=()
    y = 0.0000492*(k1o^2) + 0.0102*k1o + 0.0685
    noise = y.*randn(216)
    GFPN = GFP+noise
    k1n = mean(GFPN)
    k2n = var(GFPN)
    k3n = mean((GFPN.-k1).^3)
    k4n = mean((GFPN.-k1).^4) .- 3*k2^2
    k5n = mean((GFPN.-k1).^5) - 10*k3*k2
    slopen = ((GFPN[end]-GFPN[1])/(time[end]-time[1]))
    cvn = k2/k1
    pwn = welch_pgram(GFPN)
    psn = pwn.power[2:11]
    simn_sumstats = Float64[slopen,k1n,k2n,k3n,k4n,k5n,psn...,cvn]
    push!(simn_sumstats_list,simn_sumstats)
    push!(simn_GFP,GFPN)
end

plot!(legend=:false)
plot(time,sim_GFP,legend=:false,xlabel="Time (hrs)",ylabel="#GFP Molecules (10^6)",yscale=:log10)
plot(time,sim_GFP,legend=:false,xlabel="Time (hrs)",ylabel="#GFP Molecules (10^6)")
plot(time,simn_GFP,legend=:false,xlabel="Time (hrs)",ylabel="#GFP Molecules (10^6)")
sim_sumstats=[];
simn_sumstats=[];
sim_sumstats = reduce(vcat,transpose.(sim_sumstats_list))
simn_sumstats = reduce(vcat,transpose.(simn_sumstats_list))
params = reduce(vcat,transpose.(sim_params_list))



jldsave("0_0033_0_00033_1_25_noise_2DNABind_1.jld2"; sim_sumstats,simn_sumstats,params)

#@unpack sumstatst,paramst = jldopen("bfc.jld2")
#sim_sumstats = sumstatst
#Calculate summary statistics of the experimental data as performed on simulated data
BFC_GFP=[]
GFPex=[]
GFP=[]
ex_GFP=[]

#load and format experimental data
BFC_GFP = XLSX.readxlsx("C:/Users/Strey Lab/Documents/GitHub/Analysis_of_SGC_Simulations/StreyCats/Broken Feedback Circuit/Broken Feedback Circuit Time Traces 5_21_24 GFP.xlsx");
GFPex = BFC_GFP["Sheet1"];
GFPex = GFPex[:];
time = round.(GFPex[2:end,1].*0.33334, digits=2)
ex_GFP = GFPex[2:end,2:end]/1000000
#set up summary statics calculation loop
cols = size(ex_GFP,2);
ex_sumstats_list = [];
#plot experimental time traces
plot(ex_GFP[:,9], legend=:false,xlabel="time points",ylabel="#GFP Molecules (10^6)") #,ylims=(0,80)

#loop calculates summary statistics for experimental time traces exactly as performed on the simulated time traces
for i in 1:cols
    slope = ((ex_GFP[end,i]-ex_GFP[1,i])/(time[end]-time[1]))
    k1 = mean(ex_GFP[:,i])
    k2 = var(ex_GFP[:,i])
    k3 = mean((ex_GFP[:,i].-k1).^3)
    k4 = mean((ex_GFP[:,i].-k1).^4).-3*k2^2
    k5 = mean((ex_GFP[:,i].-k1).^5)-10*k3*k2
    pw = welch_pgram(ex_GFP[:,i])
    ps = pw.power[2:11]
    cv = k2/k1
    ex_sumstats = Float64[slope,k1,k2,k3,k4,k5,ps...,cv]
    push!(ex_sumstats_list,ex_sumstats)
end
ex_sumstats = reduce(vcat, transpose.(ex_sumstats_list));




# Visually compare the summary statistics of the simulated time traces and the experimental time traces
bins=200

# slope
min_edge = min(minimum(simn_sumstats[:,1]), minimum(ex_sumstats[:,1]))
max_edge = max(maximum(simn_sumstats[:,1]), maximum(ex_sumstats[:,1]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,1], normalize=:probability, color=:red, bins=bin_edges, alpha=0.5,label="Exp Slope")
#histogram!(sim_sumstats[:,1],bins=6000,color=:grey,alpha=0.5,label="Sim Temp Means")
histogram!(simn_sumstats[:,1], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Slope", ylabel="Probability", xlabel="Slope of GFP # (10^6)")
#savefig(MeanHist_0_18xExp_vs_Sim_Broken.png)

#Mean
min_edge = min(minimum(simn_sumstats[:,2]), minimum(ex_sumstats[:,2]))
max_edge = max(maximum(simn_sumstats[:,2]), maximum(ex_sumstats[:,2]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,2], color=:red, normalize=:probability,  bins=bin_edges, alpha=0.5, label="Exp Temp Mean")
histogram!(simn_sumstats[:,2], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Temp Mean", ylabel="Probability", xlabel="Temporal Mean of GFP # (10^6)") #,ylim=(0,0.125),xlim=(0,20)

# Variance
min_edge = min(minimum(simn_sumstats[:,3]), minimum(ex_sumstats[:,3]))
max_edge = max(maximum(simn_sumstats[:,3]), maximum(ex_sumstats[:,3]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,3], color=:red, normalize=:probability,  bins=bin_edges, alpha=0.5, label="Exp Temp Vars")
histogram!(simn_sumstats[:,3], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Temp Vars", ylabel="Probability", xlabel="Temporal Variance of GFP # (10^6)") #,ylim=(0,0.125),xlim=(0,20)

min_edge = min(minimum(simn_sumstats[:,3]), minimum(ex_sumstats[:,3]))
max_edge = max(maximum(simn_sumstats[:,3]), maximum(ex_sumstats[:,3]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,3], color=:red, normalize=:probability,  bins=bin_edges, alpha=0.5, label="Exp Temp Vars")
histogram!(simn_sumstats[:,3], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Temp Vars", ylabel="Probability", xlabel="Temporal Variance of GFP # (10^6)",ylim=(0,0.125),xlim=(0,20))

# CV
min_edge = min(minimum(simn_sumstats[:,17]), minimum(ex_sumstats[:,17]))
max_edge = max(maximum(simn_sumstats[:,17]), maximum(ex_sumstats[:,17]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,17], color=:red, normalize=:probability,  bins=200, alpha=0.5, label="Exp Temp CV")
histogram!(simn_sumstats[:,17], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Temp Vars", ylabel="Probability", xlabel="Temporal CV of GFP # (10^6)")

min_edge = min(minimum(simn_sumstats[:,17]), minimum(ex_sumstats[:,17]))
max_edge = max(maximum(simn_sumstats[:,17]), maximum(ex_sumstats[:,17]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,17], color=:red, normalize=:probability,  bins=200, alpha=0.5, label="Exp Temp CV")
histogram!(simn_sumstats[:,17], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Temp Vars", ylabel="Probability", xlabel="Temporal CV of GFP # (10^6)",ylim=(0,0.175),xlim=(0,0.5))
#k3 - skewness
min_edge = min(minimum(simn_sumstats[:,4]), minimum(ex_sumstats[:,4]))
max_edge = max(maximum(simn_sumstats[:,4]), maximum(ex_sumstats[:,4]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,4], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k3")
histogram!(simn_sumstats[:,4], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k3", ylabel="Probability", xlabel="k3") #, ylim=(0,0.20), xlim=(-30,50)

min_edge = min(minimum(simn_sumstats[:,4]), minimum(ex_sumstats[:,4]))
max_edge = max(maximum(simn_sumstats[:,4]), maximum(ex_sumstats[:,4]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,4], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k3")
histogram!(simn_sumstats[:,4], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k3", ylabel="Probability", xlabel="k3", ylim=(0,0.20), xlim=(-30,50))

#k4 - Kurtosis
min_edge = min(minimum(simn_sumstats[:,5]), minimum(ex_sumstats[:,5]))
max_edge = max(maximum(simn_sumstats[:,5]), maximum(ex_sumstats[:,5]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,5],  color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k4")
histogram!(simn_sumstats[:,5], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k4", ylabel="Probability", xlabel="k4") #, ylim=(0,0.20), xlim=(-500,250)

min_edge = min(minimum(simn_sumstats[:,5]), minimum(ex_sumstats[:,5]))
max_edge = max(maximum(simn_sumstats[:,5]), maximum(ex_sumstats[:,5]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,5],  color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k4")
histogram!(simn_sumstats[:,5], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k4", ylabel="Probability", xlabel="k4", ylim=(0,0.20), xlim=(-500,250))

#k5
min_edge = min(minimum(simn_sumstats[:,6]), minimum(ex_sumstats[:,6]))
max_edge = max(maximum(simn_sumstats[:,6]), maximum(ex_sumstats[:,6]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,6],  color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k5")
histogram!(simn_sumstats[:,6], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k5", ylabel="Probability", xlabel="k5") #, ylim=(0,0.20),xlim=(-10000,5000)

min_edge = min(minimum(simn_sumstats[:,6]), minimum(ex_sumstats[:,6]))
max_edge = max(maximum(simn_sumstats[:,6]), maximum(ex_sumstats[:,6]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,6],  color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k5")
histogram!(simn_sumstats[:,6], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k5", ylabel="Probability", xlabel="k5", ylim=(0,0.20),xlim=(-10000,5000))

#power spectrum
min_edge = min(minimum(simn_sumstats[:,7]), minimum(ex_sumstats[:,7]))
max_edge = max(maximum(simn_sumstats[:,7]), maximum(ex_sumstats[:,7]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,7], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k6")
histogram!(simn_sumstats[:,7], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k6", ylabel="Probability", xlabel="k6") #, ylim=(0,0.085),xlim=(0,10)

min_edge = min(minimum(simn_sumstats[:,7]), minimum(ex_sumstats[:,7]))
max_edge = max(maximum(simn_sumstats[:,7]), maximum(ex_sumstats[:,7]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,7], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k6")
histogram!(simn_sumstats[:,7], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k6", ylabel="Probability", xlabel="k6", ylim=(0,0.085),xlim=(0,10))


min_edge = min(minimum(simn_sumstats[:,8]), minimum(ex_sumstats[:,8]))
max_edge = max(maximum(simn_sumstats[:,8]), maximum(ex_sumstats[:,8]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,8], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k7")
histogram!(simn_sumstats[:,8], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k7", ylabel="Probability", xlabel="k7") #, ylim=(0,0.1),xlim=(0,3)

min_edge = min(minimum(simn_sumstats[:,8]), minimum(ex_sumstats[:,8]))
max_edge = max(maximum(simn_sumstats[:,8]), maximum(ex_sumstats[:,8]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,8], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k7")
histogram!(simn_sumstats[:,8], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k7", ylabel="Probability", xlabel="k7", ylim=(0,0.1),xlim=(0,3))


min_edge = min(minimum(simn_sumstats[:,9]), minimum(ex_sumstats[:,9]))
max_edge = max(maximum(simn_sumstats[:,9]), maximum(ex_sumstats[:,9]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,9], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k8")
histogram!(simn_sumstats[:,9], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k8", ylabel="Probability", xlabel="k8") #, ylim=(0,0.1),xlim=(0,2.25)

min_edge = min(minimum(simn_sumstats[:,9]), minimum(ex_sumstats[:,9]))
max_edge = max(maximum(simn_sumstats[:,9]), maximum(ex_sumstats[:,9]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,9], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k8")
histogram!(simn_sumstats[:,9], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k8", ylabel="Probability", xlabel="k8", ylim=(0,0.1),xlim=(0,2.25))


min_edge = min(minimum(simn_sumstats[:,10]), minimum(ex_sumstats[:,10]))
max_edge = max(maximum(simn_sumstats[:,10]), maximum(ex_sumstats[:,10]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,10], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k9")
histogram!(simn_sumstats[:,10], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k9", ylabel="Probability", xlabel="k9") #, ylim=(0,0.1),xlim=(0,1.75)

min_edge = min(minimum(simn_sumstats[:,10]), minimum(ex_sumstats[:,10]))
max_edge = max(maximum(simn_sumstats[:,10]), maximum(ex_sumstats[:,10]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,10], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k9")
histogram!(simn_sumstats[:,10], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k9", ylabel="Probability", xlabel="k9", ylim=(0,0.1),xlim=(0,1.75))


min_edge = min(minimum(simn_sumstats[:,11]), minimum(ex_sumstats[:,11]))
max_edge = max(maximum(simn_sumstats[:,11]), maximum(ex_sumstats[:,11]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,11], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k10")
histogram!(simn_sumstats[:,11], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k10", ylabel="Probability", xlabel="k10") #, ylim=(0,0.1),xlim=(0,1.75)

min_edge = min(minimum(simn_sumstats[:,11]), minimum(ex_sumstats[:,11]))
max_edge = max(maximum(simn_sumstats[:,11]), maximum(ex_sumstats[:,11]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,11], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k10")
histogram!(simn_sumstats[:,11], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k10", ylabel="Probability", xlabel="k10", ylim=(0,0.1),xlim=(0,1.75))


min_edge = min(minimum(simn_sumstats[:,12]), minimum(ex_sumstats[:,12]))
max_edge = max(maximum(simn_sumstats[:,12]), maximum(ex_sumstats[:,12]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,12], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k11")
histogram!(simn_sumstats[:,12], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k11", ylabel="Probability", xlabel="k11") #, ylim=(0,0.225),xlim=(0,3)

min_edge = min(minimum(simn_sumstats[:,12]), minimum(ex_sumstats[:,12]))
max_edge = max(maximum(simn_sumstats[:,12]), maximum(ex_sumstats[:,12]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,12], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k11")
histogram!(simn_sumstats[:,12], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k11", ylabel="Probability", xlabel="k11", ylim=(0,0.225),xlim=(0,3))


min_edge = min(minimum(simn_sumstats[:,13]), minimum(ex_sumstats[:,13]))
max_edge = max(maximum(simn_sumstats[:,13]), maximum(ex_sumstats[:,13]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,13], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k12")
histogram!(simn_sumstats[:,13], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k12", ylabel="Probability", xlabel="k12") #, ylim=(0,0.175),xlim=(0,1.25)

min_edge = min(minimum(simn_sumstats[:,13]), minimum(ex_sumstats[:,13]))
max_edge = max(maximum(simn_sumstats[:,13]), maximum(ex_sumstats[:,13]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,13], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k12")
histogram!(simn_sumstats[:,13], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k12", ylabel="Probability", xlabel="k12", ylim=(0,0.175),xlim=(0,1.25))


min_edge = min(minimum(simn_sumstats[:,14]), minimum(ex_sumstats[:,14]))
max_edge = max(maximum(simn_sumstats[:,14]), maximum(ex_sumstats[:,14]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,14], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k13")
histogram!(simn_sumstats[:,14], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k13", ylabel="Probability", xlabel="k13") #, ylim=(0,0.175),xlim=(0,1.25)

min_edge = min(minimum(simn_sumstats[:,14]), minimum(ex_sumstats[:,14]))
max_edge = max(maximum(simn_sumstats[:,14]), maximum(ex_sumstats[:,14]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,14], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k13")
histogram!(simn_sumstats[:,14], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k13", ylabel="Probability", xlabel="k13", ylim=(0,0.175),xlim=(0,1.25))


min_edge = min(minimum(simn_sumstats[:,15]), minimum(ex_sumstats[:,15]))
max_edge = max(maximum(simn_sumstats[:,15]), maximum(ex_sumstats[:,15]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,15], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k14")
histogram!(simn_sumstats[:,15], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k14", ylabel="Probability", xlabel="k14") #, ylim=(0,0.185),xlim=(0,1.85)

min_edge = min(minimum(simn_sumstats[:,15]), minimum(ex_sumstats[:,15]))
max_edge = max(maximum(simn_sumstats[:,15]), maximum(ex_sumstats[:,15]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,15], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k14")
histogram!(simn_sumstats[:,15], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k14", ylabel="Probability", xlabel="k14", ylim=(0,0.185),xlim=(0,1.85))


min_edge = min(minimum(simn_sumstats[:,16]), minimum(ex_sumstats[:,16]))
max_edge = max(maximum(simn_sumstats[:,16]), maximum(ex_sumstats[:,16]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,16], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label=" Exp k15")
histogram!(simn_sumstats[:,16], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k15", ylabel="Probability", xlabel="k15") #, ylim=(0,0.225),xlim=(0,2)

min_edge = min(minimum(simn_sumstats[:,16]), minimum(ex_sumstats[:,16]))
max_edge = max(maximum(simn_sumstats[:,16]), maximum(ex_sumstats[:,16]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,16], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label=" Exp k15")
histogram!(simn_sumstats[:,16], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k15", ylabel="Probability", xlabel="k15", ylim=(0,0.4),xlim=(0,2))