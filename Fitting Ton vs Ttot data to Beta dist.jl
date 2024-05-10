using Catalyst, DifferentialEquations, Plots, Interpolations
using Statistics, Distributions, StatsBase
using CSV, DataFrames, BenchmarkTools, HypothesisTests

df = CSV.read("z://Julia Programs//Promoter Statistics//Results//combined_Ton_vs_Ttot60.csv", DataFrame);
data_array = Matrix(df);

# Calculate mean and mode
mean_Ton_vs_Ttot_sims = round(mean(data_array),digits=6);
mode_Ton_vs_Ttot_sims = round(StatsBase.mode(data_array),digits=6);
stddev_Ton_vs_Ttot_sims = round(std(data_array), digits=6);

fit_beta = fit(Beta, data_array);
x_values=range(0,1, length = length(data_array));
gen_beta_dis = rand(fit_beta, length(data_array));
# Define the number of bins
num_bins = 20;

# Calculate the bin width based on the range of the data
data_range = maximum(gen_beta_dis) - minimum(gen_beta_dis);
bin_width = data_range / num_bins;

# Calculate the bin edges
bin_edges = LinRange(minimum(gen_beta_dis) - bin_width / 2, maximum(gen_beta_dis) + bin_width / 2, num_bins + 1);

# Initialize an array to store bin counts
bin_counts = zeros(Int, num_bins);

# Iterate over the data and count the number of values in each bin
for val in gen_beta_dis
    bin_index = searchsortedlast(bin_edges, val);
    if bin_index == num_bins + 1
        bin_index -= 1  # Adjust for values at the upper edge
    end
    bin_counts[bin_index] += 1
end;

# Normalize the counts to represent relative frequencies
total_observations = sum(bin_counts);
normalized_values = bin_counts / total_observations;
beta_x_values = range(0,1, length = length(normalized_values));

# Perform the Kolmogorov-Smirnov (KS) test
ks_test = ApproximateTwoSampleKSTest(gen_beta_dis, vec(data_array))
# Extract the p-value from the test result
p_value_ks = round.(pvalue(ks_test), digits=8);

Plots.histogram(data_array, xlabel="Ton / Ttot", ylabel="PDF", bins=:20, title="Fitted Histogram of Ton / Ttot\n 60 min intervals (200 sims)",
    label="Simulated Ton/Ttot Data", size = (800,800), normalize=:probability);
histogram!(gen_beta_dis,bins=:20, normalize=:probability, color=:red, alpha=0.5,label="Beta Distributed Fit");
plot!(beta_x_values, normalized_values, color=:green, label="beta fit");

vline!([mean_Ton_vs_Ttot_sims], label="Mean", color=:red, linewidth=2, linestyle=:dash);
vline!([mode_Ton_vs_Ttot_sims], label="Mode", color=:purple, linewidth=2, linestyle=:dash);


# Print mean and mode values on the histogram
annotate!([(0.83, 0.1, text("Mean: $mean_Ton_vs_Ttot_sims", :red)),
            (0.83, 0.095, text("Mode: $mode_Ton_vs_Ttot_sims", :purple)),
            (0.83, 0.090, text("Std: $stddev_Ton_vs_Ttot_sims", :blue)),
            (0.83, 0.085, text("KS Test P = $p_value_ks"))]);

plot!()


# Calculate Mean Squared Error (MSE)
mse = sum((data_array .- gen_beta_dis).^2) / length(data_array);

# Calculate Root Mean Squared Error (RMSE)
rmse = sqrt(mse)


# Display the p-value
println("Kolmogorov-Smirnov (KS) Test p-value: ", p_value_ks)

# Display the results
println("Mean Squared Error (MSE): ", mse);
println("Root Mean Squared Error (RMSE): ", rmse);
