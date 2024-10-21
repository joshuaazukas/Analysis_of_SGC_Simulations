using Catalyst, DifferentialEquations, Plots, Interpolations
using Statistics, Distributions, StatsBase
using CSV, DataFrames

@parameters kOn, kOff, kteton, ktetoff, kTr, kTl, tOn, tOff, dM, dG;
@variables t;
@species A(t), DNAoff(t), DNAon(t), DNAtetr(t), RNA(t), TetrOn(t), TetrOff(t);

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

tspan = (0.0, 600);
u0 =();
p=();
u0 = (DNAoff => 1, DNAon =>0, DNAtetr=>0, A=>1, RNA=>0, TetrOn=>0,TetrOff=>0);
p = (kOn => 0.0038*60*60, kOff => 0.0042*60*60, kteton => base_kteton*factor, ktetoff => 0.015, kTr => 375, kTl => 67.5, dM => 0.345, dG => 0.00522);

@named rn = ReactionSystem(rxs, t, [A, DNAoff, DNAon, DNAtetr, RNA, TetrOff, TetrOn], [kOn, kOff, kteton, ktetoff, kTr, kTl, tOn, tOff, dM, dG]);
rn = complete(rn)


dprob = DiscreteProblem(rn, u0, tspan, p);
jprob = JumpProblem(rn, dprob, Direct());
sol =[];
sol = solve(jprob, SSAStepper());



 t_inter = 0:0.001:600;
TetrOff_int = linear_interpolation(sol.t, sol[5,:]);
TetrOffi = TetrOff_int.(t_inter);
TetrOn_int = linear_interpolation(sol.t, sol[6,:]);

TetrOni = TetrOn_int.(t_inter)
#, xlimit=(0.0,1.0)
Tetr_sum = TetrOffi .+ TetrOni; # Doesn't work as expected
plot(sol,idxs=1)
plot(sol,idxs=2)
plot(sol,idxs=3)
plot(sol,idxs=4)
plot(sol,idxs=5)
plot(sol,idxs=6)
plot(sol,idxs=7)
#plot(t_inter, Tetr_sum)
plot(sol,idxs=7,xlimits=(800,1000))
