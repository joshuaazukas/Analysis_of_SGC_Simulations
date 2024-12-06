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
tspan = (0.0, 300);
u0 = (DNAoff => 1, DNAon => 0, A => 1, RNA => 147, GFP => 1.7098053e7); #Changed to approximately level of GFP and RNA at steady state


p = (kOn => 0.00009*60*60, kOff => 0.0001*60*60,  kTr => 310, kTl => 425, dM => 1, dG => 0.00365) #parameter set for NO FEEDBACK Circuit
ukOn = Normal(0.00009*60*60,0.2*0.00009*60*60) # standard deviation raised to 20%. estimated from Suter et al. paper that used/measured kon and koff values (no reported standard deviation)
ukOff = Normal(0.0001*60*60,0.2*0.0001*60*60)
ukTr = Normal(310,150)
# calculate gamma distribution for kTL with mean=310 and std=265
g_scale = 80^2/425  #adjusted from 265^2/500
g_shape = 425/g_scale
ukTl = Gamma(g_shape,g_scale)


# Define maximum value
max_value=[]
max_value = 575

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
GFP1 = sol[5,:][end-215:end]/1000000
#plot GFP and RNA values of the simulation
plot(sol, idxs=5, label="GFP")
plot(GFP1)
plot(sol, idxs=4,label="RNA")


#Setup for simulation loop, parameters, summary statistics, and section of GFP time trace at steady state are stored during the loop
sim_sumstats_list = Vector{Float64}[]
sim_params_list = Vector{Float64}[]
sim_GFP=Vector{Float64}[]
@time for i in 1:50000
    if (i % 100)==0
        println(i)
    end 
    kOn_new = rand(ukOn)
    kOff_new = rand(ukOff)
    kTr_new = rand(ukTr)
    kTl_new = bounded_sample(ukTl, max_value)
    u0RNA = kOn_new/(kOn_new+kOff_new)*kTr_new
    u0GFP = kTl_new*u0RNA/0.00365
    new_prob = remake(dprob; p = (kOn => kOn_new,
                            kOff => kOff_new,
                            kTr => kTr_new,
                            kTl => kTl_new),
                            u0 = (RNA => Int64(round(u0RNA)), GFP => Int64(round(u0GFP))),
                            tspan = (0.0,150.0))
    jprob = JumpProblem(rn, new_prob, Direct(); save_positions = (false, false));
    sol = solve(jprob, SSAStepper(); saveat=1/3) #changed to 0.333 hrs to match measured data
    GFP = sol[5,:][end-215:end]/1000000 #Changed to get values from indexes corresponding to past 1000hrs (steady state) now that there are 3x as many points saved in simulation output (only storing 216 data points) 
    
    k1 = mean(GFP)
    k2 = var(GFP)
    k3 = mean((GFP.-k1).^3)
    k4 = mean((GFP.-k1).^4) .- 3*k2^2
    k5 = mean((GFP.-k1).^5) - 10*k3*k2

    pw = welch_pgram(GFP)
    ps = pw.power[2:11]
    params = Float64[ kOn_new, kOff_new, kTr_new, kTl_new]
    plot!(sol, idxs=5, label="GFP $i")
    push!(sim_GFP,GFP)
    push!(sim_params_list,params)
    sim_sumstats = Float64[k1,k2,k3,k4,k5,ps...]
    push!(sim_sumstats_list,sim_sumstats)
end

plot!(legend=:false)
plot(sim_GFP,legend=:false,xlabel="time points",ylabel="#GFP Molecules (10^6)",ylims=(0,80))

sim_sumstats=[];
sim_sumstats = reduce(vcat,transpose.(sim_sumstats_list))
params = reduce(vcat,transpose.(sim_params_list))
sec_sim_traces = DataFrame([Symbol("Col$i") => vec for (i, vec) in enumerate(sim_GFP)])
#CSV.write("TestName_traces.csv", sec_sim_traces)
#jldsave("TestName.jld2"; sim_sumstats,params)