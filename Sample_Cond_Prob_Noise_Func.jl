using KernelDensity
using Distributions
using Plots
using Random
using StatsBase
using XLSX




# Estimate the joint KDE
file_path = "C://Users//Strey Lab//Documents//GitHub//Analysis_of_SGC_Simulations//StreyCats//Broken Feedback Circuit//mean_vs_noise_std.xlsx"
data = XLSX.readxlsx(file_path)
# Access the data from Sheet1
sheet = data["Sheet1"]

# Extract 'u' and 'std' columns (assuming the columns are named exactly as in your file)
u = sheet[2:end, 1]  # Get values from column 1 (u)
std = sheet[2:end, 2]  # Get values from column 2 (std)

u = convert(Array{Float64}, u)
std = convert(Array{Float64}, std)

function sample_conditional_positive_y(u::Matrix{Float64}, std::Matrix{Float64}, x_fixed::Float64)
    data = hcat(u, std)  # Combine u and std as a 2D array
    kde_result = kde(data)
    # Estimate the bivariate KDE
    kde_result = kde(data)
    
    # Extract grid points from the KDE result
    ugrid = kde_result.x
    stdgrid = kde_result.y
    density_grid = kde_result.density

    # Find the closest index in the u grid to x_fixed
    closest_u_idx = argmin(abs.(ugrid .- x_fixed))

    # Extract the density at the fixed x value for all std points
    conditional_density1 = density_grid[closest_u_idx, :]
    marginal_density_at_x = sum(conditional_density1) * (stdgrid[2] - stdgrid[1])  # Approximate integral
    # Normalize the conditional density
    conditional_density = conditional_density1 ./ marginal_density_at_x
    conditional_density /= sum(conditional_density)
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

x_fixed = 250.0        # Example fixed x value

sampled_y = sample_conditional_positive_y(u, std, x_fixed)
println("Sampled y value greater than 0: $sampled_y")
