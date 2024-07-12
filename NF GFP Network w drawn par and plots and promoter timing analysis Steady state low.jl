using Catalyst, DifferentialEquations, Plots, Random, Distributions, DataFrames, CSV, Printf, LsqFit, Interpolations
experimental_data  = CSV.read("SS Low Box Avg converted.csv", DataFrame);
observed_data = experimental_data[!,4]; #averaged GFP traces from Stead State non-induced cells, converted to approximate # of GFP molecules
observed_dataSD = experimental_data[!,5]; #SD of GFP traces, Converted by same ratio as GFP average
@parameters kOn kOff kOnt kOfft k kT deg_R deg_G;
@variables t;
@species A(t) DNA(t) A_DNA(t) DNA_T(t) A_DNA_T(t) RNA(t) GFP(t);
# Define the reaction network
rxs = [
    (@reaction kOn, A + DNA --> A_DNA), # binding of Transcription factor complex to DNA
    (@reaction kOff, A_DNA --> A + DNA), # unbinding of Transcription factor complex to DNA
    (@reaction kOnt, DNA + GFP --> DNA_T), #GFP (acting as TetR) binds to DNA, Prevents Transcription
    (@reaction kOfft, DNA_T --> DNA + GFP), # GFP unbinding DNA 
    (@reaction kOnt, A_DNA + GFP --> A_DNA_T), # GFP binding active DNA, preventing transcription
    (@reaction kOfft, A_DNA_T --> A_DNA), # GFP unbinding active DNA, allows transcription
    (@reaction k, A_DNA --> A_DNA + RNA), # Transcription of DNA to mRNA
    (@reaction kT, RNA --> RNA + GFP), # Translation of RNA to reporter protein (GFP)
    (@reaction deg_R, RNA --> 0), # mNRA Degradation
    (@reaction deg_G, GFP --> 0) # Reporter Degradation
];

tspan = (0.0, 72); # reaction time span
u0 = [A => 1, DNA => 1, A_DNA => 0, DNA_T => 0, A_DNA_T => 0, RNA => 1000, GFP => 950000];  # starting conditions

@named rn = ReactionSystem(rxs, t, [A, DNA, A_DNA, DNA_T, A_DNA_T, RNA, GFP], [kOn, kOff, kOnt, kOfft, k, kT, deg_R, deg_G]);

# Generate arrays of ~1000 values for each parameter using normal distribution
num_samples = 5000;

means = (1, 1, 1e-15, 1e12, 60, 10, 0.099, 0.03); # mean value of parameter distribution (kOn, kOff, kOnt, kOfft, k, kT, deg_R, deg_G)
std_devs = (0.1, 0.1, 1e-16, 1e11, 0.6, 0.1, 0.001, 0.001); # Standard deviation of parameter distribution
#std_devs = (0.0001, 0.0001, 0.00001, 0.0001, 0.0005, 0.00035, 0.000001, 0.000001);
# Create arrays of ~1000 values for each parameter
param_values = [rand(Normal(mean, std), num_samples) for (mean, std) in zip(means, std_devs)];
# Plot histograms for each parameter
histograms = [];
param_labels = ["kOn", "kOff", "kOnt", "kOfft", "k", "kT", "deg_R", "deg_G"];
for (i, values) in enumerate(param_values);
    hist = histogram(values, label=param_labels[i], xlabel="Parameter Value", ylabel="Frequency", bins=20);
    push!(histograms, hist);
end

# Display all histograms
plot(histograms..., layout=(2, 4), legend=:topright)
#using Subscripts, Unicode, LaTeXStrings
#khist = histogram(param_values[5], xlabel = L"k Parameter Values $(Hrs^{-1})$", ylabel = "Frequency", bins =20, color =:green, legend=false, xtickfontsize=18, ytickfontsize=18, guidefontsize=20)
#savefig(khist, "C://Users//jrazu//Desktop//Gene regulatory network Simulations//khist.png")
#degGhist = histogram(param_values[8], xlabel = L"degG Parameter Values $(Hrs^{-1})$", ylabel = "Frequency", bins =20, color =:red, legend=false, xtickfontsize=18, ytickfontsize=18, guidefontsize=20)
#savefig(degGhist, "C://Users//jrazu//Desktop//Gene regulatory network Simulations//degGhist.png")

# Arrays to store simulated GFP values and parameters for each simulation
num_simulations = 5;
threshfit = 100
saveat = 0:0.333:72;
solutions = [];
solutions_rna = [];
random_params =[];
param_sets = [];
A_DNA_sims = [];
A_DNA = [];

Random.seed!(123); # Set a seed for reproducibility
for i in 1:num_simulations;
    # Randomly choose values from the pre-generated arrays for each parameter
    random_params = [Float32(values[rand(1:num_samples)]) for values in param_values];
    push!(param_sets, random_params); # Store the parameter set for this simulation
    # Create a single DiscreteProblem and JumpProblem
    p = Dict(Symbol("kOn") => random_params[1],
             Symbol("kOff") => random_params[2],
             Symbol("kOnt") => random_params[3],
             Symbol("kOfft") => random_params[4],
             Symbol("k") => random_params[5],
             Symbol("kT") => random_params[6],
             Symbol("deg_R") => random_params[7],
             Symbol("deg_G") => random_params[8])

    dprob = DiscreteProblem(rn, u0, tspan, p)
    jprob = JumpProblem(rn, dprob, Direct())

    # Solve the problem using SSAStepper
    sol = solve(jprob, SSAStepper(), saveat=saveat)
    push!(A_DNA_sims, sol[3,:])
        # Interpolate GFP values at saveat time points
        gfp_values = []
        rna_values = []
        for t in saveat
            # Check if t is in the solution time points
            if t in sol.t
                idx = findall(x -> x == t, sol.t)[1] # Get the index of t in sol.t
                push!(gfp_values, sol[7, idx])
                push!(rna_values, sol[6, idx])
            else
                # Interpolate the GFP value at t
                interp_gfp = interpolate((sol.t,), sol[7, :], Gridded(Linear()))
                interp_rna = interpolate((sol.t,), sol[6, :], Gridded(Linear()))
                push!(gfp_values, interp_gfp(t))
                push!(rna_values, interp_rna(t))
            end
        end
    #gfp_values = sol[7,:]
    push!(solutions, gfp_values)
    push!(solutions_rna, rna_values)
end

# Reshape all GFP vectors into a single matrix with 5 columns
max_length = maximum(length, solutions)
gfp_values = hcat([vcat(vec, fill(NaN, max_length - length(vec))) for vec in solutions]...)
rna_values = hcat([vcat(vec, fill(NaN, max_length - length(vec))) for vec in solutions_rna]...)

gfp_values = gfp_values[1:end-1, :]
plot(collect(tspan[1]:0.333:tspan[2] - 0.333), gfp_values, xlabel="Time", ylabel="GFP", label="", legend=:topright);

# Plot experimental data
ps = plot(experimental_data[:, 1], observed_data, seriestype=:scatter, label="Dox Induced Experimental Time Trace", xlabel="Time (Hrs)", ylabel="# of GFP Molecules (Converted from A.u.)", color =:red, legend=:bottomright);
rna = [];
# Plot simulated data
for i in 1:5
    plot!(ps, collect(tspan[1]:0.333:tspan[2]-0.333), gfp_values[:, i], label="Example Simulation $i", linewidth=3);
end
display(ps)
plot(A_DNA_sims[1])
for i in 1:length(param_sets)
    println(param_sets[i])
end



# Calculate SSD for each simulation
ssd_values = [sum((gfp_values[:, i] .- observed_data).^2) for i in 1:num_simulations]

# Find the indices of the top 5 simulations with the minimum SSD
best_sim_indices = sortperm(ssd_values)[1:200]

# Initialize arrays to store the parameters and simulated GFP values for the top 5 simulations
top_param_sets = []
top_gfp_values = []

# Store parameters and GFP values for the top 5 simulations
for index in best_sim_indices
    # Store parameter set
    push!(top_param_sets, param_sets[index])
    
    # Store simulated GFP values
    push!(top_gfp_values, gfp_values[:, index])
end

# Plot experimental data
p = scatter(experimental_data[:, 1], observed_data, label="Dox Induced Experimental Time Trace", xlabel="Time (Hrs)", ylabel="# of GFP Molecules (Converted from A.u.)", legend=:bottomright, color = :red);

# Plot simulated data for the top 5 simulations
for (i, gfp_values) in enumerate(top_gfp_values)
    plot!(p, collect(tspan[1]:0.333:tspan[2]-0.333), gfp_values, label="Best Fit Simulation $i", linewidth=3, alpha=0.1, legend=false);
end
# Show the plot
display(p)
#l = length(top_param_sets)
# Print parameters of the top simulations
#println("Parameters of the top $l simulations:")
#for (i, params) in enumerate(top_param_sets)
#    println("Simulation $i:")
#    for (j, param) in enumerate(params)
#        println("Parameter $(param_labels[j]): $(param)")
#    end
#    println()
#end
savefig(p, "C://Users//jrazu//Desktop//Gene regulatory network Simulations//Best fit2.png")
# Define the number of bins for histograms
num_bins = 20

# Create a subplot for each parameter
histograms2 = []
for j in 1:length(param_labels)
    # Extract parameter values for the j-th parameter from the top simulations
    param_values_j = [params[j] for params in top_param_sets]
    
    # Create a histogram for the j-th parameter
    hist = histogram(param_values_j, bins=num_bins, xlabel="Parameter Value", ylabel="Frequency",color=:pink3)
    
    # Add the histogram to the list of histograms
    push!(histograms2, hist)
end

# Combine histograms into a single plot
plot(histograms2..., layout=(2, 4))

krhist = histogram(histograms2[5], xlabel = L"k Parameter Values $(Hrs^{-1})$", ylabel = "Frequency", bins =20, color =:lightgreen, legend=false, xtickfontsize=18, ytickfontsize=18, guidefontsize=20);
savefig(krhist, "C://Users//jrazu//Desktop//Gene regulatory network Simulations//krhist2.png")

degGrhist = histogram(histograms2[8], xlabel = L"degG Parameter Values $(Hrs^{-1})$", ylabel = "Frequency", bins =20, color =:lightred, legend=false, xtickfontsize=14, ytickfontsize=18, guidefontsize=20)
savefig(degGrhist, "C://Users//jrazu//Desktop//Gene regulatory network Simulations//degGrhist2.png")
