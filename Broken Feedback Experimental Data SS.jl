using Catalyst, Plots, Interpolations, JumpProcesses
using Statistics, Distributions, StatsBase
using CSV, XLSX, DataFrames, FFTW, DSP
using SymbolicIndexingInterface: parameter_values, state_values
using SciMLStructures: Tunable, replace, replace!
using JLD2, UnPack
#Retrieve Experimental Data set
BFC_GFP = XLSX.readxlsx("C:/Users/Strey Lab/Documents/GitHub/Analysis_of_SGC_Simulations/StreyCats/Broken Feedback Circuit/Broken Feedback Circuit Time Traces 5_21_24 GFP.xlsx");
GFPex = BFC_GFP["Sheet1"];
GFPex = GFPex[:];
time = round.(GFPex[2:end,1].*0.33334, digits=2)
GFP = GFPex[2:end,2:end]/1000000
GFP = GFP.*0.18
cols = size(GFP,2);
sumstats_list = [];

for i in 1:cols
    k1 = mean(GFP[:,i])
    k2 = var(GFP[:,i])
    k3 = mean((GFP[:,i].-k1).^3)
    k4 = mean((GFP[:,i].-k1).^4).-3*k2^2
    k5 = mean((GFP[:,i].-k1).^5)-10*k3*k2
    pw = welch_pgram(GFP[:,i])
    ps = pw.power[2:11]
    sumstats = Float64[k1,k2,k3,k4,k5,ps...]
    push!(sumstats_list,sumstats)
end
ex_sumstats = reduce(vcat, transpose.(sumstats_list));
jldsave("0_18xExpBroken_SumStats.jld2"; ex_sumstats)

sim_sumstats = load("sumstats_broken.jld2");
sim_sumstats = sim_sumstats["sumstatst"];

bins=200
sim_sumstatssk1 = rand(sim_sumstats[:,1],965);
min_edge = min(minimum(sim_sumstats[:,1]), minimum(ex_sumstats[:,1]))
max_edge = max(maximum(sim_sumstats[:,1]), maximum(ex_sumstats[:,1]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,1], normalize=:probability, color=:red, bins=bin_edges, alpha=0.5,label="0.18 * Exp Temp Means")
#histogram!(sim_sumstats[:,1],bins=6000,color=:grey,alpha=0.5,label="Sim Temp Means")
histogram!(sim_sumstats[:,1], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Temp Means", ylabel="Probability", xlabel="Temporal Mean of GFP # (10^6)")
#savefig(MeanHist_0_18xExp_vs_Sim_Broken.png)

sim_sumstatssk2 = rand(sim_sumstats[:,2],965);
min_edge = min(minimum(sim_sumstats[:,2]), minimum(ex_sumstats[:,2]))
max_edge = max(maximum(sim_sumstats[:,2]), maximum(ex_sumstats[:,2]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,2], color=:red, normalize=:probability,  bins=bin_edges, alpha=0.5, label="0.18 * Exp Temp Vars")
histogram!(sim_sumstats[:,2], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Temp Vars", ylabel="Probability", xlabel="Temporal Variation of GFP # (10^6)", ylim=(0,0.20))

sim_sumstatssk3 = rand(sim_sumstats[:,3],965);
min_edge = min(minimum(sim_sumstats[:,3]), minimum(ex_sumstats[:,3]))
max_edge = max(maximum(sim_sumstats[:,3]), maximum(ex_sumstats[:,3]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,3], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="0.18 * Exp k3")
histogram!(sim_sumstats[:,3], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k3", ylabel="Probability", xlabel="k3", ylim=(0,0.20))

sim_sumstatssk4 = rand(sim_sumstats[:,4],965);
min_edge = min(minimum(sim_sumstats[:,4]), minimum(ex_sumstats[:,4]))
max_edge = max(maximum(sim_sumstats[:,4]), maximum(ex_sumstats[:,4]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,4],  color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="0.18 * Exp k4")
histogram!(sim_sumstats[:,4], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k4", ylabel="Probability", xlabel="k4", ylim=(0,0.16), xlim=(-2.5,0.5))

sim_sumstatssk5 = rand(sim_sumstats[:,5],965);
min_edge = min(minimum(sim_sumstats[:,5]), minimum(ex_sumstats[:,5]))
max_edge = max(maximum(sim_sumstats[:,5]), maximum(ex_sumstats[:,5]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,5],  color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="0.18 * Exp k5")
histogram!(sim_sumstats[:,5], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k5", ylabel="Probability", xlabel="k5", ylim=(0,0.15),xlim=(-5,5))

sim_sumstatssk6 = rand(sim_sumstats[:,6],965);
min_edge = min(minimum(sim_sumstats[:,6]), minimum(ex_sumstats[:,6]))
max_edge = max(maximum(sim_sumstats[:,6]), maximum(ex_sumstats[:,6]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,6], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="0.18 * Exp k6")
histogram!(sim_sumstats[:,6], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k6", ylabel="Probability", xlabel="k6", ylim=(0,0.15),xlim=(0,0.33))

sim_sumstatssk7 = rand(sim_sumstats[:,7],965);
min_edge = min(minimum(sim_sumstats[:,7]), minimum(ex_sumstats[:,7]))
max_edge = max(maximum(sim_sumstats[:,7]), maximum(ex_sumstats[:,7]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,7], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="0.18 * Exp k7")
histogram!(sim_sumstats[:,7], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k7", ylabel="Probability", xlabel="k7", ylim=(0,0.1),xlim=(0,0.15))

sim_sumstatssk8 = rand(sim_sumstats[:,8],965);
min_edge = min(minimum(sim_sumstats[:,8]), minimum(ex_sumstats[:,8]))
max_edge = max(maximum(sim_sumstats[:,8]), maximum(ex_sumstats[:,8]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,8], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="0.18 * Exp k8")
histogram!(sim_sumstats[:,8], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k8", ylabel="Probability", xlabel="k8", ylim=(0,0.1),xlim=(0,0.15))

sim_sumstatssk9 = rand(sim_sumstats[:,9],965);
min_edge = min(minimum(sim_sumstats[:,9]), minimum(ex_sumstats[:,9]))
max_edge = max(maximum(sim_sumstats[:,9]), maximum(ex_sumstats[:,9]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,9], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="0.18 * Exp k9")
histogram!(sim_sumstats[:,9], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k9", ylabel="Probability", xlabel="k9", ylim=(0,0.1),xlim=(0,0.075))

sim_sumstatssk10 = rand(sim_sumstats[:,10],965);
min_edge = min(minimum(sim_sumstats[:,10]), minimum(ex_sumstats[:,10]))
max_edge = max(maximum(sim_sumstats[:,10]), maximum(ex_sumstats[:,10]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,10], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="0.18 * Exp k10")
histogram!(sim_sumstats[:,10], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k10", ylabel="Probability", xlabel="k10", ylim=(0,0.1),xlim=(0,0.075))

sim_sumstatssk11 = rand(sim_sumstats[:,11],965);
min_edge = min(minimum(sim_sumstats[:,11]), minimum(ex_sumstats[:,11]))
max_edge = max(maximum(sim_sumstats[:,11]), maximum(ex_sumstats[:,11]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,11], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="0.18 * Exp k11")
histogram!(sim_sumstats[:,11], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k11", ylabel="Probability", xlabel="k11", ylim=(0,0.225),xlim=(0,0.1))

sim_sumstatssk12 = rand(sim_sumstats[:,12],965);
min_edge = min(minimum(sim_sumstats[:,12]), minimum(ex_sumstats[:,12]))
max_edge = max(maximum(sim_sumstats[:,12]), maximum(ex_sumstats[:,12]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,12], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="0.18 * Exp k12")
histogram!(sim_sumstats[:,12], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k12", ylabel="Probability", xlabel="k12", ylim=(0,0.1),xlim=(0,0.075))

sim_sumstatssk13 = rand(sim_sumstats[:,13],965);
min_edge = min(minimum(sim_sumstats[:,13]), minimum(ex_sumstats[:,13]))
max_edge = max(maximum(sim_sumstats[:,13]), maximum(ex_sumstats[:,13]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,13], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="0.18 * Exp k13")
histogram!(sim_sumstats[:,13], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k13", ylabel="Probability", xlabel="k13", ylim=(0,0.175),xlim=(0,0.09))

sim_sumstatssk14 = rand(sim_sumstats[:,14],965);
min_edge = min(minimum(sim_sumstats[:,14]), minimum(ex_sumstats[:,14]))
max_edge = max(maximum(sim_sumstats[:,14]), maximum(ex_sumstats[:,14]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,14], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="0.18 * Exp k14")
histogram!(sim_sumstats[:,14], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k14", ylabel="Probability", xlabel="k14", ylim=(0,0.185),xlim=(0,0.09))

sim_sumstatssk15 = rand(sim_sumstats[:,15],965);
min_edge = min(minimum(sim_sumstats[:,15]), minimum(ex_sumstats[:,15]))
max_edge = max(maximum(sim_sumstats[:,15]), maximum(ex_sumstats[:,15]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,15], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="0.18 * Exp k15")
histogram!(sim_sumstats[:,15], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k15", ylabel="Probability", xlabel="k15", ylim=(0,0.225),xlim=(0,0.085))