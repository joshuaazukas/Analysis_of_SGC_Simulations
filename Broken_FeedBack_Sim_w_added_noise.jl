using Catalyst, Plots, Interpolations, JumpProcesses
using Statistics, Distributions, StatsBase
using CSV, XLSX, DataFrames, FFTW, DSP
using SymbolicIndexingInterface: parameter_values, state_values
using SciMLStructures: Tunable, replace, replace!
using JLD2, UnPack

@parameters kOn, kOff, kteton, ktetoff, kTr, kTl, dM, dG;
@variables t;
@species A(t), DNAoff(t), DNAon(t), RNA(t), GFP(t);

# Initial conditions and base parameters
tspan = (0.0, 2000);
u0 = (DNAoff => 1, DNAon => 0, A => 1, RNA => 0, GFP => 0); #Changed to approximately level of GFP and RNA at steady state
#base_kteton = 0.0007

p = (kOn => 0.00009*60*60, kOff => 0.0001*60*60,  kTr => 310, kTl => 425, dM => 1, dG => 0.00365) #parameter set for NO FEEDBACK Circuit
ukOn = Normal(0.00009*60*60,0.2*0.00009*60*60) # standard deviation raised to 20%. estimated from Suter et al. paper that used/measured kon and koff values (no reported standard deviation)
ukOff = Normal(0.0001*60*60,0.2*0.0001*60*60)
ukTr = Normal(310,150)
# calculate gamma distribution for kTL with mean=310 and std=265
g_scale = 80^2/425  #265^2/500
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
dist_b_kTl=[]
for i in 1:1000
    # Generate a new sample
    kTl_new = bounded_sample(ukTl, max_value)
    push!(dist_b_kTl,kTl_new)
end
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
    
# Create the ReactionSystem
@named rn = ReactionSystem(rxs, t, [A, DNAoff, DNAon, RNA, GFP], [kOn, kOff, kTr, kTl, dM, dG]);
rn = complete(rn);

# Define and solve the problem
dprob = DiscreteProblem(rn, u0, tspan, p);
jprob = JumpProblem(rn, dprob, Direct(); save_positions = (false, false));
@time sol = solve(jprob, SSAStepper(); saveat=0.333);

states = unknowns(rn)
params = parameters(rn)
GFP1 = sol[5,:][4500:4715]/1000000

plot(sol, idxs=5, label="GFP")
plot(GFP1)

k1 = mean(GFP1)

y = 0.0001*(k1^2) + 0.0068*k1 + 0.1146

noise = y.*randn(216)
histogram(noise,bins=10)

GFP1_N = GFP1+noise
plot(GFP1)
plot!(GFP1_N)