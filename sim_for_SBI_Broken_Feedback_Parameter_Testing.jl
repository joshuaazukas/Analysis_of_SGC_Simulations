using Catalyst, Plots, Interpolations, JumpProcesses
using Statistics, Distributions, StatsBase
using CSV, DataFrames, FFTW, DSP
using SymbolicIndexingInterface: parameter_values, state_values
using SciMLStructures: Tunable, replace, replace!
using JLD2, UnPack
using XLSX

@parameters kOn, kOff, kteton, ktetoff, kTr, kTl, dM, dG;
@variables t;
@species A(t), DNAoff(t), DNAon(t), RNA(t), GFP(t);

# Initial conditions and base parameters
tspan = (0.0, 2000);
u0 = (DNAoff => 1, DNAon => 0, A => 1, RNA => 0, GFP => 0); #Changed to approximately level of GFP and RNA at steady state


p = (kOn => 0.000042*60*60, kOff => 0.000052*60*60,  kTr => 310, kTl => 425, dM => 1, dG => 0.00365) #parameter set for NO FEEDBACK Circuit
ukOn = Normal(0.000042*60*60,0.2*0.000042*60*60) # standard deviation raised to 20%. estimated from Suter et al. paper that used/measured kon and koff values (no reported standard deviation)
ukOff = Normal(0.000052*60*60,0.2*0.00052*60*60)
ukTr = Normal(310,150)
# calculate gamma distribution for kTL with mean=310 and std=265
g_scale = 80^2/425  #adjusted from 265^2/500
g_shape = 425/g_scale
ukTl = Gamma(g_shape,g_scale)


# Define maximum value
max_value=[]
max_value = 750

# Function to sample with a maximum value
function bounded_sample(distribution, max_value)
    while true
        sample = rand(distribution)
        if sample <= max_value
            return sample
        end
    end
end
#create ukTl distribution with a maximum value
dist_b_kTl=[]
for i in 1:1000
    # Generate a new sample
    kTl_new = bounded_sample(ukTl, max_value)
    push!(dist_b_kTl,kTl_new)
end
#visuallize the distribution sampled from and the resultant distribution of selected kTl values, representative of what will occur during the simulation loop
kTl_x = 0:01:1000
kTl_y = pdf.(ukTl,kTl_x)
plot(kTl_x, kTl_y)
histogram(dist_b_kTl,bins=200)
# Define the reaction network
rxs = [
    (@reaction kOn, DNAoff + A --> DNAon),
    (@reaction kOff, DNAon --> DNAoff + A),
    (@reaction kTr, DNAon --> DNAon + RNA),
    (@reaction kTl, RNA --> RNA + GFP),
    (@reaction dM, RNA --> 0),
    (@reaction dG, GFP --> 0),
];
# single run before starting loop to check how the "average" trace will look based on parameters    
# Create the ReactionSystem
@named rn = ReactionSystem(rxs, t, [A, DNAoff, DNAon, RNA, GFP], [kOn, kOff, kTr, kTl, dM, dG]);
rn = complete(rn);

# Define and solve the problem
dprob = DiscreteProblem(rn, u0, tspan, p);
jprob = JumpProblem(rn, dprob, Direct(); save_positions = (false, false));
@time sol = solve(jprob, SSAStepper(); saveat=0.333);

states = unknowns(rn)
params = parameters(rn)
#define steady state section of simulation to where summary statistics will be calculated from
GFP1 = sol[5,:][4500:4715]/1000000
#plot GFP and RNA values of the simulation
plot(sol, idxs=5, label="GFP")
plot(GFP1)
plot(sol, idxs=4,label="RNA")


#Setup for simulation loop
sim_sumstats_list = Vector{Float64}[]
sim_params_list = Vector{Float64}[]
sim_GFP=Vector{Float64}[]
@time for i in 1:1000
    if (i % 100)==0
        println(i)
    end 
    kOn_new = rand(ukOn)
    kOff_new = rand(ukOff)
    kTr_new = rand(ukTr)
    kTl_new = bounded_sample(ukTl, max_value)
    new_prob = remake(dprob; p = (kOn => kOn_new,
                            kOff => kOff_new,
                            kTr => kTr_new,
                            kTl => kTl_new))
    jprob = JumpProblem(rn, new_prob, Direct(); save_positions = (false, false));
    sol = solve(jprob, SSAStepper(); saveat=0.333) #changed to 0.333 hrs to match measured data
    GFP = sol[5,:][4500:4715]/1000000 #Changed to get values from indexes corresponding to past 1000hrs (steady state) now that there are 3x as many points saved in simulation output (only storing 216 data points) 
    
    k1 = mean(GFP)
    k2 = var(GFP)
    k3 = mean((GFP.-k1).^3)
    k4 = mean((GFP.-k1).^4) .- 3*k2^2
    k5 = mean((GFP.-k1).^5) - 10*k3*k2

    pw = welch_pgram(GFPN)
    ps = pw.power[2:11]
    params = Float64[ kOn_new, kOff_new, kTr_new, kTl_new]
    plot!(sol, idxs=5, label="GFP $i")
    push!(sim_GFP,GFPN)
    push!(sim_params_list,params)
    sim_sumstats = Float64[k1,k2,k3,k4,k5,ps...]
    push!(sim_sumstats_list,sim_sumstats)
end

plot!(legend=:false)
plot(sim_GFP,legend=:false,xlabel="time points",ylabel="#GFP Molecules (10^6)",ylims=(0,80))
sim_sumstats=[];
sim_sumstats = reduce(vcat,transpose.(sim_sumstats_list))
params = reduce(vcat,transpose.(sim_params_list))

jldsave("noise_test.jld2"; sim_sumstats,params)

@unpack sim_sumstats,params = jldopen("noise_test2.jld2")

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
plot(ex_GFP, legend=:false,xlabel="time points",ylabel="#GFP Molecules (10^6)",ylims=(0,80))

#loop calculates summary statistics for experimental time traces exactly as performed on the simulated time traces
for i in 1:cols
    k1 = mean(ex_GFP[:,i])
    k2 = var(ex_GFP[:,i])
    k3 = mean((ex_GFP[:,i].-k1).^3)
    k4 = mean((ex_GFP[:,i].-k1).^4).-3*k2^2
    k5 = mean((ex_GFP[:,i].-k1).^5)-10*k3*k2
    pw = welch_pgram(ex_GFP[:,i])
    ps = pw.power[2:11]
    ex_sumstats = Float64[k1,k2,k3,k4,k5,ps...]
    push!(ex_sumstats_list,ex_sumstats)
end
ex_sumstats = reduce(vcat, transpose.(ex_sumstats_list));

# Visually compare the summary statistics of the simulated time traces and the experimental time traces
bins=200

# mean
min_edge = min(minimum(sim_sumstats[:,1]), minimum(ex_sumstats[:,1]))
max_edge = max(maximum(sim_sumstats[:,1]), maximum(ex_sumstats[:,1]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,1], normalize=:probability, color=:red, bins=bin_edges, alpha=0.5,label="Exp Temp Means")
#histogram!(sim_sumstats[:,1],bins=6000,color=:grey,alpha=0.5,label="Sim Temp Means")
histogram!(sim_sumstats[:,1], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Temp Means", ylabel="Probability", xlabel="Temporal Mean of GFP # (10^6)")
#savefig(MeanHist_0_18xExp_vs_Sim_Broken.png)

#variance
min_edge = min(minimum(sim_sumstats[:,2]), minimum(ex_sumstats[:,2]))
max_edge = max(maximum(sim_sumstats[:,2]), maximum(ex_sumstats[:,2]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,2], color=:red, normalize=:probability,  bins=bin_edges, alpha=0.5, label="Exp Temp Vars")
histogram!(sim_sumstats[:,2], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Temp Vars", ylabel="Probability", xlabel="Temporal Variation of GFP # (10^6)") #,ylim=(0,0.125),xlim=(0,20)

min_edge = min(minimum(sim_sumstats[:,2]), minimum(ex_sumstats[:,2]))
max_edge = max(maximum(sim_sumstats[:,2]), maximum(ex_sumstats[:,2]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,2], color=:red, normalize=:probability,  bins=bin_edges, alpha=0.5, label="Exp Temp Vars")
histogram!(sim_sumstats[:,2], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Temp Vars", ylabel="Probability", xlabel="Temporal Variation of GFP # (10^6)",ylim=(0,0.125),xlim=(0,20))

#k3 - skewness
min_edge = min(minimum(sim_sumstats[:,3]), minimum(ex_sumstats[:,3]))
max_edge = max(maximum(sim_sumstats[:,3]), maximum(ex_sumstats[:,3]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,3], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k3")
histogram!(sim_sumstats[:,3], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k3", ylabel="Probability", xlabel="k3") #, ylim=(0,0.20), xlim=(-30,50)

min_edge = min(minimum(sim_sumstats[:,3]), minimum(ex_sumstats[:,3]))
max_edge = max(maximum(sim_sumstats[:,3]), maximum(ex_sumstats[:,3]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,3], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k3")
histogram!(sim_sumstats[:,3], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k3", ylabel="Probability", xlabel="k3", ylim=(0,0.20), xlim=(-30,50))

#k4 - Kurtosis
min_edge = min(minimum(sim_sumstats[:,4]), minimum(ex_sumstats[:,4]))
max_edge = max(maximum(sim_sumstats[:,4]), maximum(ex_sumstats[:,4]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,4],  color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k4")
histogram!(sim_sumstats[:,4], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k4", ylabel="Probability", xlabel="k4") #, ylim=(0,0.20), xlim=(-500,250)

min_edge = min(minimum(sim_sumstats[:,4]), minimum(ex_sumstats[:,4]))
max_edge = max(maximum(sim_sumstats[:,4]), maximum(ex_sumstats[:,4]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,4],  color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k4")
histogram!(sim_sumstats[:,4], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k4", ylabel="Probability", xlabel="k4", ylim=(0,0.20), xlim=(-500,250))

#k5
min_edge = min(minimum(sim_sumstats[:,5]), minimum(ex_sumstats[:,5]))
max_edge = max(maximum(sim_sumstats[:,5]), maximum(ex_sumstats[:,5]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,5],  color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k5")
histogram!(sim_sumstats[:,5], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k5", ylabel="Probability", xlabel="k5") #, ylim=(0,0.20),xlim=(-10000,5000)

min_edge = min(minimum(sim_sumstats[:,5]), minimum(ex_sumstats[:,5]))
max_edge = max(maximum(sim_sumstats[:,5]), maximum(ex_sumstats[:,5]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,5],  color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k5")
histogram!(sim_sumstats[:,5], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k5", ylabel="Probability", xlabel="k5", ylim=(0,0.20),xlim=(-10000,5000))

#power spectrum
min_edge = min(minimum(sim_sumstats[:,6]), minimum(ex_sumstats[:,6]))
max_edge = max(maximum(sim_sumstats[:,6]), maximum(ex_sumstats[:,6]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,6], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k6")
histogram!(sim_sumstats[:,6], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k6", ylabel="Probability", xlabel="k6") #, ylim=(0,0.085),xlim=(0,10)

min_edge = min(minimum(sim_sumstats[:,6]), minimum(ex_sumstats[:,6]))
max_edge = max(maximum(sim_sumstats[:,6]), maximum(ex_sumstats[:,6]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,6], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k6")
histogram!(sim_sumstats[:,6], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k6", ylabel="Probability", xlabel="k6", ylim=(0,0.085),xlim=(0,10))


min_edge = min(minimum(sim_sumstats[:,7]), minimum(ex_sumstats[:,7]))
max_edge = max(maximum(sim_sumstats[:,7]), maximum(ex_sumstats[:,7]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,7], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k7")
histogram!(sim_sumstats[:,7], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k7", ylabel="Probability", xlabel="k7") #, ylim=(0,0.1),xlim=(0,3)

min_edge = min(minimum(sim_sumstats[:,7]), minimum(ex_sumstats[:,7]))
max_edge = max(maximum(sim_sumstats[:,7]), maximum(ex_sumstats[:,7]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,7], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k7")
histogram!(sim_sumstats[:,7], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k7", ylabel="Probability", xlabel="k7", ylim=(0,0.1),xlim=(0,3))


min_edge = min(minimum(sim_sumstats[:,8]), minimum(ex_sumstats[:,8]))
max_edge = max(maximum(sim_sumstats[:,8]), maximum(ex_sumstats[:,8]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,8], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k8")
histogram!(sim_sumstats[:,8], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k8", ylabel="Probability", xlabel="k8") #, ylim=(0,0.1),xlim=(0,2.25)

min_edge = min(minimum(sim_sumstats[:,8]), minimum(ex_sumstats[:,8]))
max_edge = max(maximum(sim_sumstats[:,8]), maximum(ex_sumstats[:,8]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,8], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k8")
histogram!(sim_sumstats[:,8], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k8", ylabel="Probability", xlabel="k8", ylim=(0,0.1),xlim=(0,2.25))


min_edge = min(minimum(sim_sumstats[:,9]), minimum(ex_sumstats[:,9]))
max_edge = max(maximum(sim_sumstats[:,9]), maximum(ex_sumstats[:,9]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,9], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k9")
histogram!(sim_sumstats[:,9], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k9", ylabel="Probability", xlabel="k9") #, ylim=(0,0.1),xlim=(0,1.75)

min_edge = min(minimum(sim_sumstats[:,9]), minimum(ex_sumstats[:,9]))
max_edge = max(maximum(sim_sumstats[:,9]), maximum(ex_sumstats[:,9]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,9], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k9")
histogram!(sim_sumstats[:,9], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k9", ylabel="Probability", xlabel="k9", ylim=(0,0.1),xlim=(0,1.75))


min_edge = min(minimum(sim_sumstats[:,10]), minimum(ex_sumstats[:,10]))
max_edge = max(maximum(sim_sumstats[:,10]), maximum(ex_sumstats[:,10]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,10], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k10")
histogram!(sim_sumstats[:,10], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k10", ylabel="Probability", xlabel="k10") #, ylim=(0,0.1),xlim=(0,1.75)

min_edge = min(minimum(sim_sumstats[:,10]), minimum(ex_sumstats[:,10]))
max_edge = max(maximum(sim_sumstats[:,10]), maximum(ex_sumstats[:,10]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,10], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k10")
histogram!(sim_sumstats[:,10], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k10", ylabel="Probability", xlabel="k10", ylim=(0,0.1),xlim=(0,1.75))


min_edge = min(minimum(sim_sumstats[:,11]), minimum(ex_sumstats[:,11]))
max_edge = max(maximum(sim_sumstats[:,11]), maximum(ex_sumstats[:,11]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,11], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k11")
histogram!(sim_sumstats[:,11], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k11", ylabel="Probability", xlabel="k11") #, ylim=(0,0.225),xlim=(0,3)

min_edge = min(minimum(sim_sumstats[:,11]), minimum(ex_sumstats[:,11]))
max_edge = max(maximum(sim_sumstats[:,11]), maximum(ex_sumstats[:,11]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,11], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k11")
histogram!(sim_sumstats[:,11], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k11", ylabel="Probability", xlabel="k11", ylim=(0,0.225),xlim=(0,3))


min_edge = min(minimum(sim_sumstats[:,12]), minimum(ex_sumstats[:,12]))
max_edge = max(maximum(sim_sumstats[:,12]), maximum(ex_sumstats[:,12]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,12], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k12")
histogram!(sim_sumstats[:,12], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k12", ylabel="Probability", xlabel="k12") #, ylim=(0,0.175),xlim=(0,1.25)

min_edge = min(minimum(sim_sumstats[:,12]), minimum(ex_sumstats[:,12]))
max_edge = max(maximum(sim_sumstats[:,12]), maximum(ex_sumstats[:,12]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,12], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k12")
histogram!(sim_sumstats[:,12], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k12", ylabel="Probability", xlabel="k12", ylim=(0,0.175),xlim=(0,1.25))


min_edge = min(minimum(sim_sumstats[:,13]), minimum(ex_sumstats[:,13]))
max_edge = max(maximum(sim_sumstats[:,13]), maximum(ex_sumstats[:,13]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,13], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k13")
histogram!(sim_sumstats[:,13], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k13", ylabel="Probability", xlabel="k13") #, ylim=(0,0.175),xlim=(0,1.25)

min_edge = min(minimum(sim_sumstats[:,13]), minimum(ex_sumstats[:,13]))
max_edge = max(maximum(sim_sumstats[:,13]), maximum(ex_sumstats[:,13]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,13], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k13")
histogram!(sim_sumstats[:,13], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k13", ylabel="Probability", xlabel="k13", ylim=(0,0.175),xlim=(0,1.25))


min_edge = min(minimum(sim_sumstats[:,14]), minimum(ex_sumstats[:,14]))
max_edge = max(maximum(sim_sumstats[:,14]), maximum(ex_sumstats[:,14]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,14], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k14")
histogram!(sim_sumstats[:,14], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k14", ylabel="Probability", xlabel="k14") #, ylim=(0,0.185),xlim=(0,1.85)

min_edge = min(minimum(sim_sumstats[:,14]), minimum(ex_sumstats[:,14]))
max_edge = max(maximum(sim_sumstats[:,14]), maximum(ex_sumstats[:,14]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,14], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k14")
histogram!(sim_sumstats[:,14], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k14", ylabel="Probability", xlabel="k14", ylim=(0,0.185),xlim=(0,1.85))


min_edge = min(minimum(sim_sumstats[:,15]), minimum(ex_sumstats[:,15]))
max_edge = max(maximum(sim_sumstats[:,15]), maximum(ex_sumstats[:,15]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,15], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label=" Exp k15")
histogram!(sim_sumstats[:,15], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k15", ylabel="Probability", xlabel="k15") #, ylim=(0,0.225),xlim=(0,2)

min_edge = min(minimum(sim_sumstats[:,15]), minimum(ex_sumstats[:,15]))
max_edge = max(maximum(sim_sumstats[:,15]), maximum(ex_sumstats[:,15]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,15], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label=" Exp k15")
histogram!(sim_sumstats[:,15], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k15", ylabel="Probability", xlabel="k15", ylim=(0,0.4),xlim=(0,2))