using Catalyst, DifferentialEquations, Plots, Interpolations
using Statistics, Distributions, StatsBase
using CSV, DataFrames

@parameters kOn, kOff, kTr, kTl, tOn, tOff, dM, dG, a;
@variables t;
@species A(t), DNAoff(t), DNAon(t), RNA(t), TetrOn(t), TetrOff(t);

rxs = [
    (@reaction kOn, DNAoff + A --> DNAon),
    (@reaction kOff, DNAon --> DNAoff + A),
    (@reaction kTr, DNAon --> DNAon + RNA),
    (@reaction kTl, RNA --> RNA + TetrOn),
    (@reaction tOn, TetrOff --> TetrOn),
    (@reaction tOff, TetrOn --> TetrOff),
    (@reaction dM, RNA --> 0),
    (@reaction dG, TetrOn --> 0),
    (@reaction dG, TetrOff --> 0)
];

tspan = (0.0, 300);
u0 =();
p=();
b = 1;
dox = 100;
u0 = (DNAoff => 1, DNAon =>0, A=>1, RNA=>0, TetrOn=>0,TetrOff=>0);
p = (kOn =>0.0042*60*60, kOff => 0.00038*60*60, kTr => 25, kTl => 12, tOn => 100, tOff => 0.01*(dox), dM => 0.12, dG => 0.01);

@named rn = ReactionSystem(rxs, t, [A, DNAoff, DNAon, RNA, TetrOff, TetrOn], [kOn, kOff, kTr, kTl, tOn, tOff, dM, dG]);


function update_parameters!(integrator)
    TetrOn_value = integrator.u[6]  
    integrator.p[kOff] = 0.00038 * 60 * 60 * (1 + integrator.p[b] * TetrOn_value)
end


function parameter_callback(integrator)
    update_parameters!(integrator)
end


cb = ContinuousCallback(parameter_callback, ())

dprob = DiscreteProblem(rn, u0, tspan, p);
jprob = JumpProblem(rn, dprob, Direct());
sol =[];
sol = solve(jprob, SSAStepper(), callback=cb);
kOff_values = [p[kOff] for (t, p) in sol]
plot!(sol.t, kOff_values, label="kOff", linestyle=:dash, linewidth=2)


t_inter = 0:0.001:300;
TetrOff_int = linear_interpolation(sol.t, sol[5,:]);
TetrOffi = TetrOff_int.(t_inter);
TetrOn_int = linear_interpolation(sol.t, sol[6,:]);
TetrOni = TetrOn_int.(t_inter)
#, xlimit=(0.0,1.0)
Tetr_sum = TetrOffi .+ TetrOni;
plot(sol,idxs=1)
plot(sol,idxs=2)
plot(sol,idxs=3)
plot(sol,idxs=4)
plot(sol,idxs=5)
plot(sol,idxs=6)
#plot(t_inter, Tetr_sum)