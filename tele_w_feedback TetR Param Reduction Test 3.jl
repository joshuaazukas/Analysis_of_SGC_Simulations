using Catalyst, DifferentialEquations, Plots, Interpolations
using Statistics, Distributions, StatsBase
using CSV, DataFrames
@parameters kOn, kOff, kteton, ktetoff, kTr, kTl, dM, dG;
@variables t;
@species A(t), DNAoff(t), DNAon(t), DNAtetr(t), RNA(t), TetrOn(t);

# Initial conditions and base parameters
tspan = (0.0, 1000);
u0 = (DNAoff => 1, DNAon => 0, DNAtetr => 0, A => 1, RNA => 0, TetrOn => 0);
base_kteton = 0.00015

# Define the different factors for kteton
kteton_factors = [1, 0.75, 0.5, 0.25, 0.1];
plot()
p=[]
# Loop over each factor, reforming the reaction network each time
for (i, factor) in enumerate(kteton_factors)
    # Adjust the value of kteton based on the current factor
    p = (kOn => 0.0042*60*60, kOff => 0.00038*60*60, kteton => base_kteton*factor, ktetoff => 0.015, kTr => 25, kTl => 12, dM => 0.08, dG => 0.01)

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
    jprob = JumpProblem(rn, dprob, Direct());
    sol = solve(jprob, SSAStepper());
    
    # Plot the result
    plot!(sol, idxs=6, label="kteton * $factor")
end

plot!()

  t_inter = 0:0.001:1000;

plot(sol,idxs=1)
plot(sol,idxs=2)
plot(sol,idxs=3)
plot(sol,idxs=4)
plot(sol,idxs=5)
plot(sol,idxs=6)
#plot(sol,idxs=7)
#plot(t_inter, Tetr_sum)
#plot(sol,idxs=7,xlimits=(800,1000))
