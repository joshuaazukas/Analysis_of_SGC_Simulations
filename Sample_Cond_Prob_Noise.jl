using KernelDensity
using Distributions
using Plots
using Random
using StatsBase
using XLSX
using CSV
using DataFrames




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

# Estimate the joint KDE
data = hcat(u, std)  # Combine u and std as a 2D array
kde_results = kde(data)
kde_result = kde(data)

# Create a grid of points to evaluate the density (for plotting)
ugrid = LinRange(minimum(sheet[2:end, 1]), maximum(sheet[2:end, 1]), 100)  # Adjust the range and grid size as needed
stdgrid = LinRange(minimum(sheet[2:end, 2]), maximum(sheet[2:end, 2]), 100)

# Create an empty array to store density values
density_values = Float64[]

# Access the density values from the KDE result
density_grid = kde_results.density

# Correct looping over the ugrid and stdgrid in Julia
for i in 1:length(ugrid)
    for j in 1:length(stdgrid)
        u_point = ugrid[i]  # Current u value in the grid
        std_point = stdgrid[j]  # Current std value in the grid
        
        # Find the closest index in the kde result grid for u_point and std_point
        u_idx = argmin(abs.(kde_result.x .- u_point))  # Find index of the closest u value
        std_idx = argmin(abs.(kde_result.y .- std_point))  # Find index of the closest std value
        
        # Retrieve the density value corresponding to this (u_point, std_point)
        push!(density_values, density_grid[u_idx, std_idx])
    end
end

density_values_matrix = reshape(density_values, length(ugrid), length(stdgrid))
# Plot the results
# Plot the results
heatmap(ugrid, stdgrid, density_values_matrix,
        xlabel="u", ylabel="std", title="Kernel Density Estimate Heatmap")


# Extract the grid points (u and std) and densities from the kde_result
ugrid = kde_result.x  # x-coordinates (u values)
stdgrid = kde_result.y  # y-coordinates (std values)
joint_density = kde_result.density  # Density values for the (u, std) grid

# Define a fixed x (u_value) for which you want to sample the conditional distribution
x_fixed = 25  # Example fixed u value

# Find the closest index in the ugrid for the x_fixed value
index_x_fixed = argmin(abs.(ugrid .- x_fixed))  # Find closest index

# Extract the joint density for the closest u value in the grid
joint_density_at_x = joint_density[index_x_fixed, :]

# Compute the marginal density for x = x_fixed by summing joint densities over all y (std)
marginal_density_at_x = sum(joint_density_at_x) * (stdgrid[2] - stdgrid[1])  # Approximate integral

# Compute the conditional density p(y | x_fixed) by dividing joint by marginal
conditional_density = joint_density_at_x ./ marginal_density_at_x

# Normalize the conditional density so it sums to 1 (turns it into a probability distribution)
conditional_density /= sum(conditional_density)

# Plot the conditional density
plot(stdgrid, conditional_density, label="Conditional Density p(y | x=$x_fixed)")

function sample_positive_y(stdgrid, conditional_density)
    while true
        sampled_y = sample(stdgrid, Weights(conditional_density))
        if sampled_y > 0
            return sampled_y
        end
    end
end

# Sample a y value based on the conditional distribution using the normalized conditional density
sampled_y = sample_positive_y(stdgrid, Weights(conditional_density))
println("Sampled y value: $sampled_y")