using Catalyst, DifferentialEquations, Plots, Interpolations
using Statistics, Distributions, StatsBase
using CSV, DataFrames

@parameters kOn, kOff, kteton, ktetoff, kTr, kTl, dM, dG;
@variables t;
@species A(t), DNAoff(t), DNAon(t), DNAtetr(t), RNA(t), TetrOn(t);

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

tspan = (0.0, 1000);
u0 =();
p=();
u0 = (DNAoff => 1, DNAon =>0, DNAtetr=>0, A=>1, RNA=>0, TetrOn=>0);
p = (kOn =>0.0042*60*60, kOff => 0.00038*60*60, kteton => 0.00015*(1),ktetoff => 0.015, kTr => 25, kTl => 12, dM => 0.08, dG => 0.01);

@named rn = ReactionSystem(rxs, t, [A, DNAoff, DNAon, DNAtetr, RNA, TetrOn], [kOn, kOff, kteton, ktetoff, kTr, kTl, dM, dG]);
rn = complete(rn)


dprob = DiscreteProblem(rn, u0, tspan, p);
jprob = JumpProblem(rn, dprob, Direct());
sol =[];
sol = solve(jprob, SSAStepper());
plot(sol,idxs=6)


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
