using OrdinaryDiffEq
using JumpProcesses
using Plots
using Distributions
using XLSX
using Polynomials
using DSP
using KernelDensity
using Statistics, StatsBase
using JLD2, UnPack

rate(u, p, t) = p.kOn * (1-u[1]) + p.kOff * u[1]
affect!(integrator) = (integrator.u[1] = (integrator.u[1] + 1 ) % 2)
crj = ConstantRateJump(rate, affect!)
dM=1.5
p = (kOn = 0.35, kOff = 0.5, kTr =300, kTl = 900, dM = dM, dG = 0.00365)
tspan = (0.0, 2000.0)

function f!(du, u, p, t)
    du[1] = 0
    du[2] = p.kTr*u[1] - p.dM*u[2]
    du[3] = p.kTl*u[2] - p.dG*u[3]
    nothing
end

uRNA = p.kOn/(p.kOn+p.kOff)*p.kTr/p.dM
uGFP = uRNA*p.kTl/p.dG
u0 = [0, uRNA, uGFP]
oprob = ODEProblem(f!, u0, tspan, p)
joprob = JumpProblem(oprob, Direct(), crj)

@time sol = solve(joprob, Tsit5(),saveat=0.333)
plot(sol[1,:], label = ["DNAA(t)"], xlabel = "t")
plot(sol[2,:], label = ["RNA"], xlabel = "t")
plot(sol[3,:][end-215:end]/1000000,label = ["GFP"], xlabel = "t")
ukOn = Normal(0.35,0.5*0.35) # standard deviation raised to 20%. estimated from Suter et al. paper that used/measured kon and koff values (no reported standard deviation)
ukOff = Normal(0.5,0.5*0.5)
#ukTrOn = Normal(0.00025*60*60,0.3*0.00025*60*60)
#ukTrOff = Normal(0.00033*60*60, 0.3*0.00033*60*60)

# calculate gamma distribution for kTL with mean=550 and std=145
g_scale = 75^2/900  #adjusted from 265^2/500  115^2/475
g_shape = 900/g_scale
ukTl = Gamma(g_shape,g_scale)

kTl_x = 0:01:1500
kTl_y = pdf.(ukTl,kTl_x)
plot(kTl_x, kTl_y, ylabel="PDF",xlabel="kTl")

kTr_g_scale = 25^2/300
kTr_g_shape =300/kTr_g_scale
ukTrg = Gamma(kTr_g_shape,kTr_g_scale)

kTr_x = 0:01:1000
kTr_y = pdf.(ukTrg ,kTr_x)
plot(kTr_x, kTr_y, label="kTrS", ylabel="PDF", xlabel="kTrS (Hr^-1)")

udG_ln = LogNormal(-5.7,0.42)

function bounded_sample_min0(distribution)
    while true
        sample = rand(distribution)
        if 0 <= sample
            return sample
        end
    end
end

# Define maximum value
min_value=250
max_value = 1250

# Function to sample with a min and maximum value
function bounded_sample(distribution, max_value, min_value)
    while true
        sample = rand(distribution)
        if min_value <= sample <= max_value
            return sample
        end
    end
end

#setup for Noise Addition
# Estimate the joint KDE
file_path = "C://Users//Strey Lab//Documents//GitHub//Analysis_of_SGC_Simulations//StreyCats//Broken Feedback Circuit//mean_vs_noise_std.xlsx" #this data is important for adding noise to simulations
#file_path = "/Users/hstrey/Documents/programming/Analysis_of_SGC_Simulations/StreyCats/Broken Feedback Circuit/mean_vs_noise_std.xlsx" #this data is important for adding noise to simulations
data = XLSX.readxlsx(file_path);
sheet = data["Sheet1"];
u = vec(sheet["A2:A966"]);
stds = vec(sheet["B2:B966"]);

u = convert(Array{Float64}, u);
stds = convert(Array{Float64}, stds);

function sample_cond_dist_std_noise(u, stds, x_fixed; cap_value=77)
    # Cap x_fixed to the maximum allowed value
    capped_x_fixed = min(x_fixed, cap_value)
    data = hcat(u, stds)
    kde_result = kde(data)
    # Create the grid
    ugrid = kde_result.x
    stdgrid = kde_result.y
    # Find the closest index to the capped x value
    u_index = argmin(abs.(ugrid .- capped_x_fixed))
    # Extract the conditional density for the capped x value
    conditional_density = kde_result.density[u_index, :]
    # Normalize the conditional density
    conditional_density /= sum(conditional_density)
    # Sample from the conditional density, ensuring the result is greater than 0
    sampled_y = sample(stdgrid, Weights(conditional_density))
    while sampled_y <= 0
        sampled_y = sample(stdgrid, Weights(conditional_density))
    end

    return sampled_y
end

#Setup for simulation loop
sim_sumstats_list = Vector{Float64}[];
simn_sumstats_list = Vector{Float64}[];
sim_params_list = Vector{Float64}[];
sim_GFP=Vector{Float64}[];
simn_GFP=Vector{Float64}[];
plot();
dir_path = "/Users/hstrey/Documents/programming/Analysis_of_SGC_Simulations/VaryAct_HHS/" #path to save data and figures
t = round.(collect(0:215).*0.33334, digits=2)
@time for i in 1:40000
    if (i % 100)==0
        println(i)
    end 
    kOn_new = bounded_sample_min0(ukOn)
    kOff_new = bounded_sample_min0(ukOff)
    kTr_new = bounded_sample(ukTrg, 1250, 100)
    kTl_new = bounded_sample(ukTl, max_value,min_value)
    dG_new = bounded_sample_min0(udG_ln)
    u0RNA_new = kOn_new/(kOn_new+kOff_new)*kTr_new
    u0GFP_new = kTl_new*u0RNA_new/dG_new
    u0_new = [0,u0RNA_new,u0GFP_new]
    p1 = (kOn = kOn_new,kOff = kOff_new,kTr = kTr_new,kTl = kTl_new,dG = dG_new, dM = dM)
    #new_prob = remake(oprob; p = (kOn = kOn_new,
                            #kOff = kOff_new,
                            #kTr = kTr_new,
                            #kTl = kTl_new,
                            #dG = dG_new),
                            #u0 = (0,u0RNA_new,u0GFP_new),
                            #tspan = (0.0,500.0))
    #oprob = ODEProblem(f!, u0, tspan, p)                     
    #new_prob = remake(oprob; p= p1)
                            
    tspan = (0.0,2000.0)
    oprob = ODEProblem(f!, u0_new, tspan, p1)
    joprob = JumpProblem(oprob, Direct(), crj)
    sol = solve(joprob, Tsit5(),saveat=0.333) #changed to 0.333 hrs to match measured data
    GFP = sol[3,:][end-215:end]/1000000 #Changed to get values from indexes corresponding to past 1000hrs (steady state) now that there are 3x as many points saved in simulation output (only storing 216 data points) 
    # summary statistics of output w/o noise 
    slope = Polynomials.fit(t, GFP, 1)[1]
    k1 = mean(GFP)
    k2 = var(GFP)
    k3 = mean((GFP.-k1).^3)
    k4 = mean((GFP.-k1).^4) .- 3*k2^2
    k5 = mean((GFP.-k1).^5) - 10*k3*k2
    cv = k2/k1
    pw = welch_pgram(GFP)
    ps = pw.power[2:11]
    params = Float64[ kOn_new, kOff_new, kTr_new, kTl_new, dG_new]
    push!(sim_GFP,GFP)
    push!(sim_params_list,params)
    sim_sumstats = Float64[slope,k1,k2,k3,k4,k5,ps...]
    push!(sim_sumstats_list,sim_sumstats)

    k1o =()
    k1o = mean(GFP)

    #add noise
    y=()
    y = sample_cond_dist_std_noise(u, stds, k1o)
    noise = y.*randn(216)
    GFPN = GFP+noise
    k1n = mean(GFPN)
    k2n = var(GFPN)
    k3n = mean((GFPN.-k1n).^3)
    k4n = mean((GFPN.-k1n).^4) .- 3*k2n^2
    k5n = mean((GFPN.-k1n).^5) - 10*k3n*k2n
    slopen = Polynomials.fit(t, GFPN, 1)[1]
    cvn = k2n/k1n
    pwn = welch_pgram(GFPN)
    psn = pwn.power[2:11]
    simn_sumstats = Float64[slopen,k1n,k2n,k3n,k4n,k5n,psn...]
    plot!(sol, idxs=3, label="GFP $i")
    push!(simn_sumstats_list,simn_sumstats)
    push!(simn_GFP,GFPN)
end

# plot and save simulation outputs
plot!(legend=:false)
savefig(joinpath(dir_path, "FullSimOutput.png"))
plot(t,sim_GFP,legend=:false,xlabel="Time (hrs)",ylabel="#GFP Molecules (10^6)")
savefig(joinpath(dir_path, "SimOutput.png"))
plot(t,simn_GFP,legend=:false,xlabel="Time (hrs)",ylabel="#GFP Molecules (10^6)")
savefig(joinpath(dir_path, "SimOutputN.png"))

sim_sumstats=[];
simn_sumstats=[];
sim_sumstats = reduce(vcat,transpose.(sim_sumstats_list))
simn_sumstats = reduce(vcat,transpose.(simn_sumstats_list))
params = reduce(vcat,transpose.(sim_params_list));

histogram(params[:,1], bins=200,xlabel="kOn",ylabel="counts",label="Mean:$(round(mean(params[:,1]), digits=2))")
#savefig(joinpath(dir_path, "kOn.png"))
histogram(params[:,2], bins=200,xlabel="kOff",ylabel="counts",label="Mean:$(round(mean(params[:,2]), digits=2))")
#savefig(joinpath(dir_path, "kOff.png"))
histogram(params[:,3], bins=200,xlabel="kTr",ylabel="counts",label="Mean:$(round(mean(params[:,3]), digits=2))")
histogram(params[:,4], bins=200,xlabel="kTl",ylabel="counts",label="Mean:$(round(mean(params[:,4]), digits=2))")
histogram(params[:,5], bins=200,xlabel="dG",ylabel="counts",label="Mean:$(round(mean(params[:,5]), digits=4))")

@unpack ex_sumstats = jldopen("BrokenFeedbackExperimental.jld2");
# compare summary statistic distributions of experimental and simulated + noise data
plot(ex_sumstats[:,3],ex_sumstats[:,2], seriestype=:scatter, color=:red, xlabel="Variance", ylabel="Mean",label="exp")
plot!(simn_sumstats[:,3],simn_sumstats[:,2], seriestype=:scatter, color=:grey, xlabel="Variance", ylabel="Mean",label="sim")

plot(ex_sumstats[:,3],ex_sumstats[:,1], seriestype=:scatter, color=:red, xlabel="Variance", ylabel="Slope",label="exp")
plot!(simn_sumstats[:,3],simn_sumstats[:,1], seriestype=:scatter, color=:grey, xlabel="Variance", ylabel="Slope",label="sim")

bins=200
# slope
min_edge = min(minimum(simn_sumstats[:,1]), minimum(ex_sumstats[:,1]))
max_edge = max(maximum(simn_sumstats[:,1]), maximum(ex_sumstats[:,1]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,1], normalize=:probability, color=:red, bins=bin_edges, alpha=0.5,label="Exp Slope")
#histogram!(sim_sumstats[:,1],bins=6000,color=:grey,alpha=0.5,label="Sim Temp Means")
histogram!(simn_sumstats[:,1], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Slope", ylabel="Probability", xlabel="Slope of GFP # (10^6)")

#Mean
min_edge = min(minimum(simn_sumstats[:,2]), minimum(ex_sumstats[:,2]))
max_edge = max(maximum(simn_sumstats[:,2]), maximum(ex_sumstats[:,2]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,2], color=:red, normalize=:probability,  bins=bin_edges, alpha=0.5, label="Exp Temp Mean")
histogram!(simn_sumstats[:,2], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Temp Mean", ylabel="Probability", xlabel="Temporal Mean of GFP # (10^6)")

# Variance
min_edge = min(minimum(simn_sumstats[:,3]), minimum(ex_sumstats[:,3]))
max_edge = max(maximum(simn_sumstats[:,3]), maximum(ex_sumstats[:,3]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,3], color=:red, normalize=:probability,  bins=bin_edges, alpha=0.5, label="Exp Temp Vars")
histogram!(simn_sumstats[:,3], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Temp Vars", ylabel="Probability", xlabel="Temporal Variance of GFP # (10^6)") #,ylim=(0,0.125),xlim=(0,20)

min_edge = min(minimum(simn_sumstats[:,3]), minimum(ex_sumstats[:,3]))
max_edge = max(maximum(simn_sumstats[:,3]), maximum(ex_sumstats[:,3]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,3], color=:red, normalize=:probability,  bins=bin_edges, alpha=0.5, label="Exp Temp Vars")
histogram!(simn_sumstats[:,3], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim Temp Vars", ylabel="Probability", xlabel="Temporal Variance of GFP # (10^6)",ylim=(0,0.125),xlim=(0,20))

#k3 - skewness
min_edge = min(minimum(simn_sumstats[:,4]), minimum(ex_sumstats[:,4]))
max_edge = max(maximum(simn_sumstats[:,4]), maximum(ex_sumstats[:,4]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,4], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k3")
histogram!(simn_sumstats[:,4], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k3", ylabel="Probability", xlabel="k3") #, ylim=(0,0.20), xlim=(-30,50)
savefig(joinpath(dir_path, "SkewExp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,4]), minimum(ex_sumstats[:,4]))
max_edge = max(maximum(simn_sumstats[:,4]), maximum(ex_sumstats[:,4]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,4], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k3")
histogram!(simn_sumstats[:,4], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k3", ylabel="Probability", xlabel="k3", ylim=(0,0.20), xlim=(-30,50))
savefig(joinpath(dir_path, "SkewzExp_vs_Sim_Broken.png"))
#k4 - Kurtosis
min_edge = min(minimum(simn_sumstats[:,5]), minimum(ex_sumstats[:,5]))
max_edge = max(maximum(simn_sumstats[:,5]), maximum(ex_sumstats[:,5]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,5],  color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k4")
histogram!(simn_sumstats[:,5], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k4", ylabel="Probability", xlabel="k4") #, ylim=(0,0.20), xlim=(-500,250)
savefig(joinpath(dir_path, "KurtExp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,5]), minimum(ex_sumstats[:,5]))
max_edge = max(maximum(simn_sumstats[:,5]), maximum(ex_sumstats[:,5]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,5],  color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k4")
histogram!(simn_sumstats[:,5], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k4", ylabel="Probability", xlabel="k4", ylim=(0,0.20), xlim=(-500,250))
savefig(joinpath(dir_path, "KurtzExp_vs_Sim_Broken.png"))
#k5
min_edge = min(minimum(simn_sumstats[:,6]), minimum(ex_sumstats[:,6]))
max_edge = max(maximum(simn_sumstats[:,6]), maximum(ex_sumstats[:,6]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,6],  color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k5")
histogram!(simn_sumstats[:,6], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k5", ylabel="Probability", xlabel="k5") #, ylim=(0,0.20),xlim=(-10000,5000)
savefig(joinpath(dir_path, "K5Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,6]), minimum(ex_sumstats[:,6]))
max_edge = max(maximum(simn_sumstats[:,6]), maximum(ex_sumstats[:,6]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,6],  color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k5")
histogram!(simn_sumstats[:,6], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k5", ylabel="Probability", xlabel="k5", ylim=(0,0.20),xlim=(-10000,5000))
savefig(joinpath(dir_path, "K5zExp_vs_Sim_Broken.png"))
#power spectrum
min_edge = min(minimum(simn_sumstats[:,7]), minimum(ex_sumstats[:,7]))
max_edge = max(maximum(simn_sumstats[:,7]), maximum(ex_sumstats[:,7]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,7], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k6")
histogram!(simn_sumstats[:,7], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k6", ylabel="Probability", xlabel="k6") #, ylim=(0,0.085),xlim=(0,10)
savefig(joinpath(dir_path, "PS1Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,7]), minimum(ex_sumstats[:,7]))
max_edge = max(maximum(simn_sumstats[:,7]), maximum(ex_sumstats[:,7]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,7], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k6")
histogram!(simn_sumstats[:,7], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k6", ylabel="Probability", xlabel="k6", ylim=(0,0.09),xlim=(0,15))
savefig(joinpath(dir_path, "PS1zExp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,8]), minimum(ex_sumstats[:,8]))
max_edge = max(maximum(simn_sumstats[:,8]), maximum(ex_sumstats[:,8]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,8], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k7")
histogram!(simn_sumstats[:,8], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k7", ylabel="Probability", xlabel="k7") #, ylim=(0,0.1),xlim=(0,3)
savefig(joinpath(dir_path, "PS2Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,8]), minimum(ex_sumstats[:,8]))
max_edge = max(maximum(simn_sumstats[:,8]), maximum(ex_sumstats[:,8]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,8], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k7")
histogram!(simn_sumstats[:,8], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k7", ylabel="Probability", xlabel="k7", ylim=(0,0.11),xlim=(0,6))
savefig(joinpath(dir_path, "PS2zExp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,9]), minimum(ex_sumstats[:,9]))
max_edge = max(maximum(simn_sumstats[:,9]), maximum(ex_sumstats[:,9]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,9], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k8")
histogram!(simn_sumstats[:,9], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k8", ylabel="Probability", xlabel="k8") #, ylim=(0,0.1),xlim=(0,2.25)
savefig(joinpath(dir_path, "PS3Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,9]), minimum(ex_sumstats[:,9]))
max_edge = max(maximum(simn_sumstats[:,9]), maximum(ex_sumstats[:,9]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,9], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k8")
histogram!(simn_sumstats[:,9], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k8", ylabel="Probability", xlabel="k8", ylim=(0,0.5),xlim=(0,8))
savefig(joinpath(dir_path, "PS3zExp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,10]), minimum(ex_sumstats[:,10]))
max_edge = max(maximum(simn_sumstats[:,10]), maximum(ex_sumstats[:,10]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,10], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k9")
histogram!(simn_sumstats[:,10], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k9", ylabel="Probability", xlabel="k9") #, ylim=(0,0.1),xlim=(0,1.75)
savefig(joinpath(dir_path, "PS4Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,10]), minimum(ex_sumstats[:,10]))
max_edge = max(maximum(simn_sumstats[:,10]), maximum(ex_sumstats[:,10]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,10], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k9")
histogram!(simn_sumstats[:,10], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k9", ylabel="Probability", xlabel="k9", ylim=(0,0.8),xlim=(0,8))
savefig(joinpath(dir_path, "PS4zExp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,11]), minimum(ex_sumstats[:,11]))
max_edge = max(maximum(simn_sumstats[:,11]), maximum(ex_sumstats[:,11]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,11], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k10")
histogram!(simn_sumstats[:,11], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k10", ylabel="Probability", xlabel="k10") #, ylim=(0,0.1),xlim=(0,1.75)
savefig(joinpath(dir_path, "PS5Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,11]), minimum(ex_sumstats[:,11]))
max_edge = max(maximum(simn_sumstats[:,11]), maximum(ex_sumstats[:,11]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,11], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k10")
histogram!(simn_sumstats[:,11], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k10", ylabel="Probability", xlabel="k10", ylim=(0,0.85),xlim=(0,8))
savefig(joinpath(dir_path, "PS5zExp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,12]), minimum(ex_sumstats[:,12]))
max_edge = max(maximum(simn_sumstats[:,12]), maximum(ex_sumstats[:,12]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,12], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k11")
histogram!(simn_sumstats[:,12], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k11", ylabel="Probability", xlabel="k11") #, ylim=(0,0.225),xlim=(0,3)
savefig(joinpath(dir_path, "PS6Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,12]), minimum(ex_sumstats[:,12]))
max_edge = max(maximum(simn_sumstats[:,12]), maximum(ex_sumstats[:,12]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,12], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k11")
histogram!(simn_sumstats[:,12], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k11", ylabel="Probability", xlabel="k11", ylim=(0,0.8),xlim=(0,8))
savefig(joinpath(dir_path, "PS6zExp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,13]), minimum(ex_sumstats[:,13]))
max_edge = max(maximum(simn_sumstats[:,13]), maximum(ex_sumstats[:,13]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,13], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k12")
histogram!(simn_sumstats[:,13], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k12", ylabel="Probability", xlabel="k12") #, ylim=(0,0.175),xlim=(0,1.25)
savefig(joinpath(dir_path, "PS7Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,13]), minimum(ex_sumstats[:,13]))
max_edge = max(maximum(simn_sumstats[:,13]), maximum(ex_sumstats[:,13]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,13], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k12")
histogram!(simn_sumstats[:,13], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k12", ylabel="Probability", xlabel="k12", ylim=(0,0.9),xlim=(0,8))
savefig(joinpath(dir_path, "PS7zExp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,14]), minimum(ex_sumstats[:,14]))
max_edge = max(maximum(simn_sumstats[:,14]), maximum(ex_sumstats[:,14]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,14], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k13")
histogram!(simn_sumstats[:,14], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k13", ylabel="Probability", xlabel="k13") #, ylim=(0,0.175),xlim=(0,1.25)
savefig(joinpath(dir_path, "PS8Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,14]), minimum(ex_sumstats[:,14]))
max_edge = max(maximum(simn_sumstats[:,14]), maximum(ex_sumstats[:,14]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,14], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k13")
histogram!(simn_sumstats[:,14], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k13", ylabel="Probability", xlabel="k13", ylim=(0,0.8),xlim=(0,8))
savefig(joinpath(dir_path, "PS8zExp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,15]), minimum(ex_sumstats[:,15]))
max_edge = max(maximum(simn_sumstats[:,15]), maximum(ex_sumstats[:,15]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,15], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k14")
histogram!(simn_sumstats[:,15], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k14", ylabel="Probability", xlabel="k14") #, ylim=(0,0.185),xlim=(0,1.85)
savefig(joinpath(dir_path, "PS9Exp_vs_Sim_Broken.png"))
min_edge = min(minimum(simn_sumstats[:,15]), minimum(ex_sumstats[:,15]))
max_edge = max(maximum(simn_sumstats[:,15]), maximum(ex_sumstats[:,15]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,15], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label="Exp k14")
histogram!(simn_sumstats[:,15], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k14", ylabel="Probability", xlabel="k14", ylim=(0,0.8),xlim=(0,8))
savefig(joinpath(dir_path, "PS9zExp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,16]), minimum(ex_sumstats[:,16]))
max_edge = max(maximum(simn_sumstats[:,16]), maximum(ex_sumstats[:,16]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,16], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label=" Exp k15")
histogram!(simn_sumstats[:,16], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k15", ylabel="Probability", xlabel="k15") #, ylim=(0,0.225),xlim=(0,2)
savefig(joinpath(dir_path, "PS10Exp_vs_Sim_Broken.png"))

min_edge = min(minimum(simn_sumstats[:,16]), minimum(ex_sumstats[:,16]))
max_edge = max(maximum(simn_sumstats[:,16]), maximum(ex_sumstats[:,16]))
bin_edges = range(min_edge, max_edge, length=bins)
histogram(ex_sumstats[:,16], color=:red, normalize=:probability, bins=bin_edges, alpha=0.5, label=" Exp k15")
histogram!(simn_sumstats[:,16], color=:grey, normalize=:probability, bins=bin_edges, alpha=0.5, label="Sim k15", ylabel="Probability", xlabel="k15", ylim=(0,0.8),xlim=(0,8))
savefig(joinpath(dir_path, "PS10zExp_vs_Sim_Broken.png"))