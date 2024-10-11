using Catalyst, DifferentialEquations, Plots, Interpolations
using Statistics, Distributions, StatsBase
using CSV, DataFrames, FFTW, DSP
using SymbolicIndexingInterface: parameter_values, state_values
using SciMLStructures: Tunable, replace, replace!
@parameters kOn, kOff, kteton, ktetoff, kTr, kTl, dM, dG;
@variables t;
@species A(t), DNAoff(t), DNAon(t), DNAtetr(t), RNA(t), TetrOn(t);

# Initial conditions and base parameters
tspan = (0.0, 1124);
u0 = (DNAoff => 1, DNAon => 0, DNAtetr => 0, A => 1, RNA => 0, TetrOn => 0);
base_kteton = 0.0007

# Define the different factors for kteton
kteton_factors=[1]
#kteton_factors = [1, 0.75, 0.5, 0.25, 0.1];

p = (kOn => 0.00042*60*60, kOff => 0.0042*60*60, kteton => 0.0001016273663, ktetoff => 54, kTr => 375, kTl => 67.5, dM => 0.345, dG => 0.1)
ukOn = Normal(0.00042*60*60,0.1*0.00042*60*60)
ukOff = Normal(0.0042*60*60,0.1*0.0042*60*60)
ukteton = Normal(0.0001016273663,0.1*0.0001016273663)
uktetoff = Normal(54,0.1*54)
ukTr = Normal(375,0.1*375)
ukTl = Normal(67.5,0.1*67.5)


# Define the reaction network
rxs = [
    (@reaction kOn, DNAoff + A --> DNAon),
    (@reaction kteton, DNAoff + TetrOn --> DNAtetr),
    (@reaction ktetoff, DNAtetr --> DNAoff + TetrOn),
    (@reaction kOff, DNAon --> DNAoff + A),
    (@reaction kTr, DNAon --> DNAon + RNA),
    (@reaction kTl, RNA --> RNA + TetrOn),
    (@reaction dM, RNA --> 0),
    (@reaction dG, TetrOn --> 0),
];
    
# Create the ReactionSystem
@named rn = ReactionSystem(rxs, t, [A, DNAoff, DNAon, DNAtetr, RNA, TetrOn], [kOn, kOff, kteton, ktetoff, kTr, kTl, dM, dG]);
rn = complete(rn);

# Define and solve the problem
dprob = DiscreteProblem(rn, u0, tspan, p);
jprob = JumpProblem(rn, dprob, Direct(); save_positions = (false, false));
@time sol = solve(jprob, SSAStepper(); saveat=1.0);

states = unknowns(rn)
params = parameters(rn)

plot(sol, idxs=6, label="kteton",xlims=(200,1124))
savefig("newpara.png")

tetr = sol[6,:][102:1125]/10000

k1 = mean(tetr)
k2 = var(tetr)
k3 = mean((tetr.-k1).^3)
k4 = mean((tetr.-k1).^4) .- 3*k2^2
k5 = mean((tetr.-k1).^5) - 10*k3*k2

tetr_fft = rfft(tetr)
tetr_auto = abs.((tetr_fft .* conj.(tetr_fft)))

plot(tetr_auto[2:end])

pw = welch_pgram(tetr)
plot(pw.power[2:end])

sumstats_list = Vector{Float64}[]
params_list = Vector{Float64}[]
@time for i in 1:10
    kOn_new = rand(ukOn)
    kOff_new = rand(ukOff)
    kteton_new = rand(ukteton)
    ktetoff_new = rand(uktetoff)
    kTr_new = rand(ukTr)
    kTl_new = rand(ukTl)
    new_prob = remake(dprob; p = (kOn => kOn_new,
                            kOff => kOff_new,
                            kteton => kteton_new,
                            ktetoff => ktetoff_new,
                            kTr => kTr_new,
                            kTl => kTl_new))
    jprob = JumpProblem(rn, dprob, Direct(); save_positions = (false, false));
    sol = solve(jprob, SSAStepper(); saveat=1.0)
    tetr = sol[6,:][102:1125]/10000

    k1 = mean(tetr)
    k2 = var(tetr)
    k3 = mean((tetr.-k1).^3)
    k4 = mean((tetr.-k1).^4) .- 3*k2^2
    k5 = mean((tetr.-k1).^5) - 10*k3*k2

    pw = welch_pgram(tetr)
    ps = pw.power[2:11]
    params = Float64[ kOn_new, kOff_new, kteton_new, ktetoff_new, kTr_new, kTl_new]
    push!(params_list,params)
    sumstats = Float64[k1,k2,k3,k4,k5,ps...]
    push!(sumstats_list,sumstats)
end

sumstats = reduce(vcat,transpose.(sumstats_list))
params = reduce(vcat,transpose.(params_list))

