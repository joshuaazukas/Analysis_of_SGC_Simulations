using Catalyst, Plots, Interpolations, JumpProcesses
using Statistics, Distributions, StatsBase
using CSV, DataFrames, FFTW, DSP
using SymbolicIndexingInterface: parameter_values, state_values
using SciMLStructures: Tunable, replace, replace!
using JLD2, UnPack

@parameters kOn, kOff, kteton, ktetoff, kTr, kTl, dM, dG;
@variables t;
@species A(t), DNAoff(t), DNAon(t), RNA(t), GFP(t);

# Initial conditions and base parameters
tspan = (0.0, 500);
u0 = (DNAoff => 1, DNAon => 0, A => 1, RNA => 35, GFP => 2800000); #Changed to approximately level of GFP and RNA at steady state
#base_kteton = 0.0007

p = (kOn => 0.00042*60*60, kOff => 0.0042*60*60,  kTr => 375, kTl => 310, dM => 1, dG => 0.00365) #parameter set for NO FEEDBACK Circuit
ukOn = Normal(0.00042*60*60,0.2*0.00042*60*60) # standard deviation raised to 20%. estimated from Suter et al. paper that used/measured kon and koff values (no reported standard deviation)
ukOff = Normal(0.0042*60*60,0.2*0.0042*60*60)
ukTr = Normal(375,125)
# calculate gamma distribution for kTL with mean=310 and std=265
g_scale = 265^2/310
g_shape = 310/g_scale
ukTl = Gamma(g_shape,g_scale)

kTl_x = 0:01:600
kTl_y = pdf.(ukTl,kTl_x)
plot(kTl_x, kTl_y)

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

plot(sol, idxs=5, label="GFP")
#plot(sol, idxs=4, label="RNA",xlims=(0,50))

plot(sol, idxs=3, label="DNAon",xlims=(0,2000))
savefig("newpara.png")

GFP = sol[5,:][1288:1503]/1000000 #Changed to get values from indexes corresponding to past 1000hrs (steady state) now that there are 3x as many points saved in simulation output
plot(GFP) #look at stored GFP trace over this time period
k1 = mean(GFP)
k2 = var(GFP)
k3 = mean((GFP.-k1).^3)
k4 = mean((GFP.-k1).^4) .- 3*k2^2
k5 = mean((GFP.-k1).^5) - 10*k3*k2

GFP_fft = rfft(GFP)
GFP_auto = abs.((GFP_fft .* conj.(GFP_fft)))

plot(GFP_auto[2:end])

pw = welch_pgram(GFP)
plot(pw.power[2:end])

sumstats_list = Vector{Float64}[]
params_list = Vector{Float64}[]
@time for i in 1:50000
    if (i % 100)==0
        println(i)
    end 
    kOn_new = rand(ukOn)
    kOff_new = rand(ukOff)
    kTr_new = rand(ukTr)
    kTl_new = rand(ukTl)
    new_prob = remake(dprob; p = (kOn => kOn_new,
                            kOff => kOff_new,
                            kTr => kTr_new,
                            kTl => kTl_new))
    jprob = JumpProblem(rn, new_prob, Direct(); save_positions = (false, false));
    sol = solve(jprob, SSAStepper(); saveat=0.333) #changed to 0.333 hrs to match measured data
    GFP = sol[5,:][1288:1503]/1000000 #Changed to get values from indexes corresponding to past 1000hrs (steady state) now that there are 3x as many points saved in simulation output (only storing 216 data points) 
    #RNA = sol[4,:][1288:1503]

    k1 = mean(GFP)
    k2 = var(GFP)
    k3 = mean((GFP.-k1).^3)
    k4 = mean((GFP.-k1).^4) .- 3*k2^2
    k5 = mean((GFP.-k1).^5) - 10*k3*k2

    #println(k1," ",mean(RNA))
    pw = welch_pgram(GFP)
    ps = pw.power[2:11]
    params = Float64[ kOn_new, kOff_new, kTr_new, kTl_new]
    #plot!(sol, idxs=5, label="GFP $i")
    push!(params_list,params)
    sumstats = Float64[k1,k2,k3,k4,k5,ps...]
    push!(sumstats_list,sumstats)
end
#plot!()
sumstats = reduce(vcat,transpose.(sumstats_list))
params = reduce(vcat,transpose.(params_list))

jldsave("4thbroken.jld2"; sumstats,params)

@unpack sumstats,params = jldopen("4thbroken.jld2")
sumstatst = sumstats
paramst = params

@unpack sumstats,params = jldopen("1stbroken.jld2")
sumstatst = vcat(sumstatst, sumstats)
paramst = vcat(paramst,params)

@unpack sumstats,params = jldopen("2ndbroken.jld2")
sumstatst = vcat(sumstatst, sumstats)
paramst = vcat(paramst,params)

@unpack sumstats,params = jldopen("3rdbroken.jld2")
sumstatst = vcat(sumstatst, sumstats)
paramst = vcat(paramst,params)

file = jldopen("sumstats_broken.jld2","w")
@pack! file = sumstatst,paramst
close(file)
f = load("sumstats_broken.jld2")
keys(f)