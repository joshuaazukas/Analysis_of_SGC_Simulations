using Catalyst, DifferentialEquations, Plots, Interpolations
using Statistics, Distributions, StatsBase
using CSV, DataFrames, BenchmarkTools, Base.Threads

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
u0 = [A => 1, DNA => 1, A_DNA => 0, DNA_T => 0, A_DNA_T => 0, RNA => 100000, GFP => 1000000];  # starting conditions
p = [kOn => 10, kOff => 10, kOnt=> 10, kOfft=> 10, k => 3, kT => 1.5, deg_R => 0.03, deg_G => 0.008]

@named rn = ReactionSystem(rxs, t, [A, DNA, A_DNA, DNA_T, A_DNA_T, RNA, GFP], [kOn, kOff, kOnt, kOfft, k, kT, deg_R, deg_G]);

num_simulations = 5;
T_of_sims = [];
A_DNA_sims = [];
GFP_sims = [];
A_DNA = [];
T_A_DNA = [];
on_state_durations = [];
on_state_durations_sims = [];
total_on_state_duration = [];
total_on_state_durations = [];
total_on_state_durations_sims = [];
sub_matrices = [];
sub_matrix =[];
T_A_DNA =[];
total_on_state_durations_sims = Vector{Vector{Float64}}(undef, num_simulations)

# Preallocate memory for sub_matrices and on_state_durations
sub_matrices = Vector{Tuple{Vector{Float64}, Vector{Int}}}(undef, num_simulations)
on_state_durations = Vector{Float64}(undef, num_simulations)
bin_width = 0.3333;
@time begin
    for i in 1:num_simulations;
        GFP=[];
        T_A_DNA =[];
        A_DNA = [];
        dprob = DiscreteProblem(rn, u0, tspan, p)
        jprob = JumpProblem(rn, dprob, Direct())

        sol = solve(jprob, SSAStepper())
        GFP = sol[7,:]
        A_DNA = sol[3,:]
        T_A_DNA = (sol.t, A_DNA);
        push!(GFP_sims, GFP)
        # Initialize an empty array to store the individual matrices
        sub_matrices = [];
        on_state_durations = [];
        total_on_state_durations =[];
        indices = [];
        bin_start =[];
        bin_end =[];
        time_points =[];
        states=[];
        sub_matrix=[];
        # Define the bin width
        bin_width = bin_width;
        # Iterate over each bin
        for bin_start in 0:bin_width:(72)
            # Define the end of the current bin
            bin_end = bin_start + bin_width
            
            # Find indices corresponding to the current bin
            
            indices = findall(x -> bin_start <= x < bin_end, T_A_DNA[1])
            
            # Extract the corresponding time points and states
            time_points = T_A_DNA[1][indices]
            states = T_A_DNA[2][indices]
            
            # Create a submatrix for the current bin
            sub_matrix = (time_points, states)
            
            # Store the submatrix in the array
            push!(sub_matrices, sub_matrix)
        end
        print(length(sub_matrices))

        # Initialize an array to store durations of the 1 state between transitions from 0 to 1 and 1 to 0
        on_state_durations = [];
        first_transition = [];
        last_transition = [];
        transitions_1_to_0 = [];
        transitions_0_to_1 = [];
        time_points =[];
        states=[];
        total_on_state_duration = [];

        #start loop here
        for x in 1:length(sub_matrices)
            on_state_durations = [];
            first_transition = [];
            last_transition = [];
            transitions_1_to_0 = [];
            transitions_0_to_1 = [];
            time_points = [];
            states = [];

            # Extract the time points and states arrays from the first submatrix
            time_points = sub_matrices[x][1];
            states = sub_matrices[x][2];
            # Find indices where the state transitions occur
            transition_indices_1_to_0 = findall(diff(states) .== -1);
            transition_indices_0_to_1 = findall(diff(states) .== 1);
            for idx in transition_indices_1_to_0
                push!(transitions_1_to_0, time_points[idx])
            end;
            for idx in transition_indices_0_to_1
                push!(transitions_0_to_1, time_points[idx])
            end;
            
            # Calculate total on state duration for the current submatrix
            start_time = time_points[1]
            end_time = time_points[end]

            #if sub_matrix starts w promoter on and ends w promoter off        
            if states[1] == 1 && states[end] == 0 #if sub_matrix starts w promoter on and ends w promoter off
                first_transition = transitions_1_to_0[1] - start_time
                push!(on_state_durations, first_transition)
                        
                for i in 1:length(transitions_0_to_1)
                    duration = (transitions_1_to_0[i + 1] - transitions_0_to_1[i])
                    push!(on_state_durations, duration)
                end
                
            #if sub_matrix starts w promoter off and ends w promoter on
            elseif states[end] == 1 && states[1] == 0 
                for i in 1:length(transitions_1_to_0)
                    duration = transitions_1_to_0[i] - transitions_0_to_1[i]
                    push!(on_state_durations, duration)
                end
                        
                last_transition = end_time - transitions_0_to_1[end]
                push!(on_state_durations, last_transition)
            
            #if promoter is on at start and end of submatrix, but has transitions 
            elseif states[1] == 1 && states[end] == 1 && length(transition_indices_0_to_1) > 0 
                first_transition = transitions_1_to_0[1] - start_time
                push!(on_state_durations, first_transition)
                        
                for i in 1:length(transitions_0_to_1) - 1
                    duration = (transitions_1_to_0[i + 1] - transitions_0_to_1[i])
                    push!(on_state_durations, duration)
                end
                        
                last_transition = end_time - transitions_0_to_1[end]
                push!(on_state_durations, last_transition)
            
            #if Promoter is off at start and end of sub_matrix, but has transitions
            elseif states[1] == 0 && states[end] == 0 && length(transition_indices_0_to_1) > 0
                for i in 1:length(transitions_0_to_1)
                    duration = transitions_1_to_0[i] - transitions_0_to_1[i]
                    push!(on_state_durations, duration)
                end
            
            #if promoter is on for the full durations of sub_matrix
            elseif length(transition_indices_0_to_1) == 0 && length(transition_indices_1_to_0) == 0 && states[1] == 1 && states[end] == 1
                duration = end_time - start_time
                push!(on_state_durations, duration)
                    
            #if promoter is off for the duration of the sub_matrix
            elseif length(transition_indices_0_to_1) == 0 && length(transition_indices_1_to_0) == 0 && states[1] == 0 && states[end] == 0
                duration = 0
                push!(on_state_durations, duration)
            
            end
                    
            total_on_state_duration = round.(sum(on_state_durations), digits = 6)
            push!(total_on_state_durations, total_on_state_duration)
            
        end
        push!(total_on_state_durations_sims, total_on_state_durations)
    end
end
Ton_vs_Ttot_per_simulation = []

# Loop over each array in the matrix
for durations_array in total_on_state_durations_sims
    # Calculate Ton_vs_Ttot for the current simulation
    current_Ton_vs_Ttot = round.((durations_array / bin_width), digits=6)
    
    # Append the result to the array
    push!(Ton_vs_Ttot_per_simulation, current_Ton_vs_Ttot)
end
combined_Ton_vs_Ttot = []

# Loop over each array in Ton_vs_Ttot_per_simulation
for Ton_vs_Ttot_array in Ton_vs_Ttot_per_simulation
    # Concatenate the current array with the combined vector
    combined_Ton_vs_Ttot = vcat(combined_Ton_vs_Ttot, Ton_vs_Ttot_array)
end
# Calculate mean and mode
mean_Ton_vs_Ttot_sims = round(mean(combined_Ton_vs_Ttot),digits=6);
mode_Ton_vs_Ttot_sims = round(StatsBase.mode(combined_Ton_vs_Ttot),digits=6);
stddev_Ton_vs_Ttot_sims = round(std(combined_Ton_vs_Ttot), digits=6);


# Plot histogram
histogram(combined_Ton_vs_Ttot, xlabel="Ton / Ttot", ylabel="Frequency", bins=:20, title="Histogram of Ton / Ttot\n 20 min intervals (50 sims)",
    label=false, size = (800,800), legend=false);
vline!([mean_Ton_vs_Ttot_sims], label="Mean", color=:red, linewidth=2, linestyle=:dash);
vline!([mode_Ton_vs_Ttot_sims], label="Mode", color=:purple, linewidth=2, linestyle=:dash);

# Print mean and mode values on the histogram
annotate!([(0.9, 750, text("Mean: $mean_Ton_vs_Ttot_sims", :red)),
            (0.9, 725, text("Mode: $mode_Ton_vs_Ttot_sims", :purple)),
            (0.9, 700, text("Std: $stddev_Ton_vs_Ttot_sims", :blue))]);

plot!()
#CSV.write("C://Users//Strey Lab//Desktop//Github Reps//single-cell-circuits-main//Promoter Statistics Analysis//A_DNA_SS High 1.csv", DataFrame(df_A_DNA))

# Initialize arrays to store results for each simulation
total_on_times = Float64[]
fractions_on_time = Float64[]

# Loop over each simulation
for i in 1:num_simulations
    # Initialize variables to track state and time
    in_1_state = false
    start_time = NaN
    on_durations_minutes = Float64[]

    # Iterate over each time point and value in the A_DNA array
    for (j, val) in enumerate(A_DNA_sims[i])
        # Check if the value is 1 and if we're not already in the 1 state
        if val == 1 && !in_1_state
            in_1_state = true
            start_time = T_of_sims[i][j]  # Record the start time in simulation time
        # Check if the value is 0 and if we're currently in the 1 state
        elseif val == 0 && in_1_state
            # Calculate the duration in 1 state and store it in minutes
            duration_minutes = (T_of_sims[i][j] - start_time) * 60  # Convert duration to minutes
            push!(on_durations_minutes, duration_minutes)
            in_1_state = false
        end
    end

    # If the simulation ends in the 1 state, calculate the duration until the end in minutes
    if in_1_state
        duration_minutes = (tspan[2] - start_time) * 60  # Convert duration to minutes
        push!(on_durations_minutes, duration_minutes)
    end

    # Calculate total time the promoter is in the "On" state in minutes for this simulation
    total_on_time_minutes = sum(on_durations_minutes)
    push!(total_on_times, total_on_time_minutes)

    # Calculate total simulation time in minutes for this simulation
    total_simulation_time_minutes = (tspan[2] - tspan[1]) * 60  # Convert hours to minutes

    # Calculate the fraction of time the promoter is in the "On" state compared to the total simulation time for this simulation
    fraction_on_time = total_on_time_minutes / total_simulation_time_minutes
    push!(fractions_on_time, fraction_on_time)
end

# Print the results
println("Total On Times for Each Simulation: ", total_on_times)
println("Fractions of Time Promoter is On for Each Simulation: ", fractions_on_time)

# Initialize arrays to store results for each simulation
average_reactions_per_intervals = Float64[]
stddev_reactions_per_intervals = Float64[]

# Loop over each simulation
for i in 1:num_simulations
    # Define the bin width
    bin_width = 0.333
    # Initialize an array to store the counts for each bin
    bin_counts = []

    # Define the range of values
    min_value = minimum(T_of_sims[i])
    max_value = maximum(T_of_sims[i])

    # Iterate over each bin
    for bin_start in min_value:bin_width:max_value
        # Define the end of the current bin
        bin_end = min(bin_start + bin_width, max_value)
        
        # Count the number of points within the current bin
        count = sum((T_of_sims[i] .>= bin_start) .& (T_of_sims[i] .< bin_end))
        
        # Store the count for the current bin
        push!(bin_counts, count)
    end

    # Compute the average number of reactions per interval for this simulation
    average_reactions_per_interval = mean(bin_counts)
    push!(average_reactions_per_intervals, average_reactions_per_interval)

    # Compute the standard deviation of the number of reactions per interval for this simulation
    stddev_reactions_per_interval = std(bin_counts)
    push!(stddev_reactions_per_intervals, stddev_reactions_per_interval)
end

# Print the results
println("Average Reactions per Interval for Each Simulation: ", average_reactions_per_intervals)
println("Standard Deviation of Reactions per Interval for Each Simulation: ", stddev_reactions_per_intervals)

T_A_DNA_sims = [];
for y in 1:(num_simulations)
    push!(T_A_DNA_sims, A_DNA_sims[y])
    push!(T_A_DNA_sims, T_of_sims[y])
end
# Initialize an array to store the pairs of data matrices
comb_T_A_DNA = []

# Iterate over each simulation
for i in 1:num_simulations
    # Create a matrix for the current simulation pair
    pair_matrix = [A_DNA_sims[i] T_of_sims[i]]
        
    # Store the pair matrix along with its label
    push!(comb_T_A_DNA, pair_matrix)
end

# Initialize an empty array to store the individual matrices
sub_matrices = [];

# Define the bin width
bin_width = 1;

# Loop over each pair of data matrices
for pair_matrix in comb_T_A_DNA
    # Initialize an array to store the submatrices for the current pair
    pair_sub_matrices = []
    
    # Iterate over each bin
    for bin_start in 0:bin_width:72
        # Define the end of the current bin
        bin_end = bin_start + bin_width
        
        # Find indices corresponding to the current bin
        indices = findall(x -> bin_start <= x < bin_end, pair_matrix[:, 2])
        
        # Extract the corresponding time points and states
        time_points = pair_matrix[indices, 2]
        states = pair_matrix[indices, 1]
        
        # Create a submatrix for the current bin
        sub_matrix = hcat(time_points, states)
        
        # Store the submatrix in the array for the current pair
        push!(pair_sub_matrices, sub_matrix)
    end
    
    # Store the submatrices for the current pair in the overall array
    push!(sub_matrices, pair_sub_matrices)
end

# Initialize an array to store durations of the 1 state between transitions from 0 to 1 and 1 to 0 for each pair
total_on_state_durations = [];

# Iterate over each pair of submatrices
for pair_sub_matrices in sub_matrices
    # Initialize an array to store the total on state durations for the current pair
    pair_total_on_state_durations = []
    
    # Iterate over each submatrix for the current pair
    for sub_matrix in pair_sub_matrices
        # Initialize variables for tracking state and time
        on_state_durations = []
        first_transition = []
        last_transition = []
        transitions_1_to_0 = []
        transitions_0_to_1 = []
        time_points = []
        states = []
        
        # Extract the time points and states arrays from the submatrix
        time_points = sub_matrix[:, 1]
        states = sub_matrix[:, 2]
        
        # Find indices where the state transitions occur
        transition_indices_1_to_0 = findall(diff(states) .== -1)
        transition_indices_0_to_1 = findall(diff(states) .== 1)
        
        for idx in transition_indices_1_to_0
            push!(transitions_1_to_0, time_points[idx])
        end
        
        for idx in transition_indices_0_to_1
            push!(transitions_0_to_1, time_points[idx])
        end
        
        # Calculate total on state duration for the current submatrix
        start_time = time_points[1]
        end_time = time_points[end]
        
        if states[1] == 1 && states[end] == 0
            first_transition = transitions_1_to_0[1] - start_time
            push!(on_state_durations, first_transition)
            
            for i in 1:length(transitions_0_to_1)
                duration = (transitions_1_to_0[i + 1] - transitions_0_to_1[i])
                push!(on_state_durations, duration)
            end
        elseif states[end] == 1 && states[1] == 0
            for i in 1:length(transitions_1_to_0)
                duration = transitions_1_to_0[i] - transitions_0_to_1[i]
                push!(on_state_durations, duration)
            end
            
            last_transition = end_time - transitions_0_to_1[end]
            push!(on_state_durations, last_transition)
        elseif states[1] == 1 && states[end] == 1 && length(transition_indices_0_to_1) > 0
            first_transition = transitions_1_to_0[1] - start_time
            push!(on_state_durations, first_transition)
            
            for i in 1:length(transitions_0_to_1) - 1
                duration = (transitions_1_to_0[i + 1] - transitions_0_to_1[i])
                push!(on_state_durations, duration)
            end
            
            last_transition = end_time - transitions_0_to_1[end]
            push!(on_state_durations, last_transition)
        elseif states[1] == 0 && states[end] == 0 && length(transition_indices_0_to_1) > 0
            for i in 1:length(transitions_0_to_1)
                duration = transitions_1_to_0[i] - transitions_0_to_1[i]
                push!(on_state_durations, duration)
            end
        elseif length(transition_indices_0_to_1) == 0 && length(transition_indices_1_to_0) == 0 && states[1] == 1 && states[end] == 1
            duration = end_time - start_time
            push!(on_state_durations, duration)
        elseif length(transition_indices_0_to_1) == 0 && length(transition_indices_1_to_0) == 0 && states[1] == 0 && states[end] == 0
            duration = 0
            push!(on_state_durations, duration)
        end
        
        total_on_state_duration = sum(on_state_durations)
        push!(pair_total_on_state_durations, total_on_state_duration)
    end
    
    # Store the total on state durations for the current pair
    push!(total_on_state_durations, pair_total_on_state_durations)
end
# Initialize an array to store rounded total on state durations for each simulation
rounded_total_on_state_durations = []

# Loop over each array in the matrix
for durations_array in total_on_state_durations
    # Round the durations and append them to the new array
    push!(rounded_total_on_state_durations, round.(durations_array, digits=6))
end

# Initialize an array to store Ton_vs_Ttot for each simulation
Ton_vs_Ttot_per_simulation = []

# Loop over each array in the matrix
for durations_array in rounded_total_on_state_durations
    # Calculate Ton_vs_Ttot for the current simulation
    current_Ton_vs_Ttot = round.((durations_array / bin_width), digits=6)
    
    # Append the result to the array
    push!(Ton_vs_Ttot_per_simulation, current_Ton_vs_Ttot)
end
# Initialize an empty vector to store combined Ton_vs_Ttot values
combined_Ton_vs_Ttot = []

# Loop over each array in Ton_vs_Ttot_per_simulation
for Ton_vs_Ttot_array in Ton_vs_Ttot_per_simulation
    # Concatenate the current array with the combined vector
    combined_Ton_vs_Ttot = vcat(combined_Ton_vs_Ttot, Ton_vs_Ttot_array)
end

# Calculate mean and mode
mean_Ton_vs_Ttot_sims = round(mean(combined_Ton_vs_Ttot),digits=6);
mode_Ton_vs_Ttot_sims = round(StatsBase.mode(combined_Ton_vs_Ttot),digits=6);
stddev_Ton_vs_Ttot_sims = round(std(combined_Ton_vs_Ttot), digits=6);


# Plot histogram
histogram(combined_Ton_vs_Ttot, xlabel="Ton / Ttot", ylabel="Frequency", bins=:20, title="Histogram of Ton / Ttot\n 60 min intervals (50 sims)",
    label=false, size = (800,800), legend=false);
vline!([mean_Ton_vs_Ttot_sims], label="Mean", color=:red, linewidth=2, linestyle=:dash);
vline!([mode_Ton_vs_Ttot_sims], label="Mode", color=:purple, linewidth=2, linestyle=:dash);

# Print mean and mode values on the histogram
annotate!([(0.9, 3000, text("Mean: $mean_Ton_vs_Ttot_sims", :red)),
            (0.9, 2500, text("Mode: $mode_Ton_vs_Ttot_sims", :purple)),
            (0.9, 2000, text("Std: $stddev_Ton_vs_Ttot_sims", :blue))]);

plot!()
maximum(combined_Ton_vs_Ttot)