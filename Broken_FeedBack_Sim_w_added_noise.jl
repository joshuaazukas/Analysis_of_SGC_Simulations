using Catalyst, Plots, Interpolations, JumpProcesses
using Statistics, Distributions, StatsBase
using CSV, XLSX, DataFrames, FFTW, DSP
using SymbolicIndexingInterface: parameter_values, state_values
using SciMLStructures: Tunable, replace, replace!
using JLD2, UnPack

@parameters kOn, kOff, kteton, ktetoff, kTr, kTl, dM, dG;
@variables t;
@species A(t), DNAoff(t), DNAon(t), RNA(t), GFP(t);

# Initial conditions and parameters
tspan = (0.0, 2000);
u0 = (DNAoff => 1, DNAon => 0, A => 1, RNA => 0, GFP => 0); #Changed to approximately level of GFP and RNA at steady state

p = (kOn => 0.000042*60*60, kOff => 0.000052*60*60,  kTr => 310, kTl => 425, dM => 1, dG => 0.00365) #parameter set for NO FEEDBACK Circuit
ukOn = Normal(0.000042*60*60,0.2*0.000042*60*60) # standard deviation raised to 20%. estimated from Suter et al. paper that used/measured kon and koff values (no reported standard deviation)
ukOff = Normal(0.000052*60*60,0.2*0.00052*60*60)
ukTr = Normal(310,150)
# calculate gamma distribution for kTL with mean=310 and std=265
g_scale = 80^2/425  #265^2/500
g_shape = 425/g_scale
ukTl = Gamma(g_shape,g_scale)


# Define maximum value cutoff
max_value=[]
max_value = 750

# Function to sample from a distribution with a maximum value cutoff
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

plot(sol, idxs=5, label="GFP")
plot(GFP1)
#calculate temporal mean of simulated trace
k1 = mean(GFP1)
#calculate std of noise to apply to simulated traces
y = 0.0001*(k1^2) + 0.0068*k1 + 0.1146
#apply noise to the simulated trace
noise = y.*randn(216)
GFP1_N = GFP1+noise
#compare simulated trace to simulated trace with added noise
plot(GFP1)
plot!(GFP1_N)

#Setup for simulation loop
sim_sumstats_list = Vector{Float64}[]
sim_params_list = Vector{Float64}[]
sim_GFP=Vector{Float64}[]
sim_GFPN=Vector{Float64}[]
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
    #RNA = sol[4,:][1288:1503]
    push!(sim_GFP,GFP)
    k1o =()
    k1o = mean(GFP)

    #add noise
    y=()
    y = 0.0000492*(k1o^2) + 0.0102*k1o + 0.0685
    noise = y.*randn(216)
    GFPN = GFP+noise
    k1 = mean(GFPN)
    k2 = var(GFPN)
    k3 = mean((GFPN.-k1).^3)
    k4 = mean((GFPN.-k1).^4) .- 3*k2^2
    k5 = mean((GFPN.-k1).^5) - 10*k3*k2

    #println(k1," ",mean(RNA))
    pw = welch_pgram(GFPN)
    ps = pw.power[2:11]
    params = Float64[ kOn_new, kOff_new, kTr_new, kTl_new]
    plot!(sol, idxs=5, label="GFP $i")
    push!(sim_GFPN,GFPN)
    push!(sim_params_list,params)
    sim_sumstats = Float64[k1,k2,k3,k4,k5,ps...]
    push!(sim_sumstats_list,sim_sumstats)
end

 plot!(legend=:false)

plot(sim_GFP,legend=:false,xlabel="time points",ylabel="#GFP Molecules (10^6)",ylims=(0,80))
plot(sim_GFPN,legend=:false,xlabel="time points",ylabel="#GFP Molecules (10^6)",ylims=(0,80))
sim_sumstats=[];
sim_sumstats = reduce(vcat,transpose.(sim_sumstats_list))
params = reduce(vcat,transpose.(sim_params_list))

jldsave("noise_test2.jld2"; sim_sumstats,params)