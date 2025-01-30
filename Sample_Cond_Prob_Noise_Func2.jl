using KernelDensity
using Distributions
using Plots
using Random
using StatsBase
using XLSX
using Statistics




# Estimate the joint KDE
file_path = "C://Users//Strey Lab//Documents//GitHub//Analysis_of_SGC_Simulations//StreyCats//Broken Feedback Circuit//mean_vs_noise_std.xlsx"
data1 = XLSX.readxlsx(file_path)
# Access the data from Sheet1
sheet = data1["Sheet1"]
data=[]
# Extract 'u' and 'std' columns (assuming the columns are named exactly as in your file)
u = sheet[2:end, 1]  # Get values from column 1 (u)
std = sheet[2:end, 2]  # Get values from column 2 (std)

u = convert(Array{Float64}, u)
std = convert(Array{Float64}, std)

function sample_conditional_with_extrapolation_and_plot(u::Matrix{Float64}, std::Matrix{Float64}, x_fixed::Float64)
    data = hcat(u, std)  # Combine u and std as a 2D array
    kde_result = kde(data)
    # Extract grid points from the KDE result
    ugrid = kde_result.x
    stdgrid = kde_result.y
    density_grid = kde_result.density

    # Check if x_fixed is within the range of ugrid
    if x_fixed < minimum(ugrid) || x_fixed > maximum(ugrid)
        println("x_fixed = $x_fixed is outside the range of ugrid (min: $(minimum(ugrid)), max: $(maximum(ugrid))).")

        # Fit a Gaussian to the tail of the data
        mean_std = mean(vec(std))  # Ensure std is a 1D vector
        std_std = std(vec(std))    # Compute standard deviation of the vector
        tail_std = vec(std)[vec(std) .> mean_std + std_std]  # Extract right tail values  # Using right tail
        if isempty(tail_std)
            error("No tail data available for extrapolation.")
        end
        tail_dist = fit(Normal, tail_std)
        
        # Plot the Gaussian distribution
        x_vals = range(mean(tail_dist) - 4 * std(tail_dist), mean(tail_dist) + 4 * std(tail_dist), length=500)
        y_vals = pdf.(tail_dist, x_vals)
        plot(x_vals, y_vals, label="Fitted Gaussian", title="Extrapolated Tail Distribution", xlabel="std", ylabel="Density")
        
        # Sample from the Gaussian until a positive value is obtained
        while true
            sampled_y = rand(tail_dist)
            if sampled_y > 0
                return sampled_y
            end
        end
    else
        # Find the closest index in the u grid to x_fixed
        closest_u_idx = argmin(abs.(ugrid .- x_fixed))

        # Extract the density at the fixed x value for all std points
        conditional_density = density_grid[closest_u_idx, :]

        # Normalize the conditional density
        conditional_density = conditional_density ./ sum(conditional_density)

        # Plot the conditional density
        plot(stdgrid, conditional_density, label="Conditional Density at x = $x_fixed", 
             title="Conditional Density Distribution", xlabel="std", ylabel="Density")
        
        # Function to sample a positive y value
        function sample_positive_y(stdgrid, conditional_density)
            while true
                sampled_y = sample(stdgrid, Weights(conditional_density))
                if sampled_y > 0
                    return sampled_y
                end
            end
        end

        # Sample a positive y value
        sampled_y_positive = sample_positive_y(stdgrid, conditional_density)
        return sampled_y_positive
    end
end


x_fixed = 250.0        # Example fixed x value

sampled_y = sample_conditional_with_extrapolation_and_plot(u, std, x_fixed)
println("Sampled y value greater than 0: $sampled_y")
