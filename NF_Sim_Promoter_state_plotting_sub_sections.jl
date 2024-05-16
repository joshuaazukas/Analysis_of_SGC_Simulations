using Catalyst, DifferentialEquations, Plots, Interpolations
using Statistics, Distributions, StatsBase
using CSV, DataFrames

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

dprob = DiscreteProblem(rn, u0, tspan, p)
jprob = JumpProblem(rn, dprob, Direct())

sol = solve(jprob, SSAStepper())

A_DNA = sol[3,:];

# Initialize an empty array to store durations in 1 state in minutes
on_durations_minutes = Float64[];

# Initialize variables to track state and time
in_1_state = false;
start_time = NaN;

# Iterate over each time point and value in the A_DNA array
for (i, val) in enumerate(A_DNA)
    # Check if the value is 1 and if we're not already in the 1 state
    if val == 1 && !in_1_state
        in_1_state = true
        start_time = sol.t[i]  # Record the start time in simulation time
    # Check if the value is 0 and if we're currently in the 1 state
    elseif val == 0 && in_1_state
        # Calculate the duration in 1 state and store it in minutes
        duration_minutes = (sol.t[i] - start_time) * 60  # Convert duration to minutes
        push!(on_durations_minutes, duration_minutes)
        in_1_state = false
    end
end

# If the simulation ends in the 1 state, calculate the duration until the end in minutes
if in_1_state
    duration_minutes = (tspan[2] - start_time) * 60  # Convert duration to minutes
    push!(on_durations_minutes, duration_minutes)
end

# Plot a histogram of on durations in minutes
histogram(on_durations_minutes, xlabel="Duration Promoter is On (minutes)", ylabel="Frequency", bins=20,
 legend=false, title = "Pormoter On State Durations for entire 72 hours", xticks = 12, xlim = (0,60), ylim = (0,240))

# Calculate total time the promoter is in the "On" state in minutes
total_on_time_minutes = sum(on_durations_minutes)
# Calculate total simulation time in minutes
total_simulation_time_minutes = (tspan[2] - tspan[1]) * 60  # Convert hours to minutes
# Calculate the fraction of time the promoter is in the "On" state compared to the total simulation time
fraction_on_time = total_on_time_minutes / total_simulation_time_minutes
# Print the fraction
println("Fraction of time promoter is in On state: ", fraction_on_time)
# Define the bin width
bin_width = 0.333;
# Initialize an array to store the counts for each bin
bin_counts = [];
# Define the range of values
min_value = minimum(sol.t);
max_value = maximum(sol.t);
# Iterate over each bin
for bin_start in min_value:bin_width:max_value
    # Define the end of the current bin
    bin_end = min(bin_start + bin_width, max_value)
    
    # Count the number of points within the current bin
    count = sum((sol.t .>= bin_start) .& (sol.t .< bin_end))
    
    # Store the count for the current bin
    push!(bin_counts, count)
end;
# Define the bin centers for plotting
bin_centers = collect(min_value : bin_width : max_value);
# Plot the histogram
bar(bin_centers, bin_counts, xlabel="Time [hrs]", ylabel="Counts", title="number of points per 20 min of 72 hour sim",
 legend=false,)
# Compute the average number of reactions per interval
average_reactions_per_interval = mean(bin_counts);

plot(sol.t[1:50], title = "time of simulation", xlabel="# of Events", ylabel="Time [Hrs]")
# Compute the standard deviation of the number of reactions per interval
stddev_reactions_per_interval = std(bin_counts);
# Display the result as the average plus or minus the standard deviation
println("Average ± Standard deviation: ", average_reactions_per_interval, " ± ", stddev_reactions_per_interval)

T_A_DNA = (sol.t,A_DNA);
plot(sol)
plot(sol.t,A_DNA, title="Promoter State for Entire 72hr simulation",
 xlabel = "time [hrs]", legend=false, ylabel="Promoter State (1 = active)")
# Initialize an empty array to store the individual matrices
sub_matrices = [];

# Define the bin width
bin_width = 0.33333;
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

int_in_min = round.(bin_width*60, digits=2)
subm_plot = 2
length(sub_matrices[1][1])
plot(sub_matrices[subm_plot], title= "Sub-Inteval$subm_plot\n $int_in_min Min Window of 72hr simulation", 
xlabel = "time [hrs]", legend=false, ylabel="Promoter State (1 = active)")

plot((sub_matrices[subm_plot][1]*60),sub_matrices[subm_plot][2], title= "Sub-Inteval$subm_plot\n $int_in_min Min Window of 72hr simulation", 
xlabel = "time [min]", legend=false, ylabel="Promoter State (1 = active)")

# Initialize an array to store durations of the 1 state between transitions from 0 to 1 and 1 to 0
total_on_state_durations = [];
on_state_durations = [];
first_transition = [];
last_transition = [];
transitions_1_to_0 = [];
transitions_0_to_1 = [];
time_points =[]
states=[]
#start loop here
for x in 1:length(sub_matrices)
    on_state_durations = []
    first_transition = [];
    last_transition = [];
    transitions_1_to_0 = [];
    transitions_0_to_1 = [];
    time_points =[]
    states=[]
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
    println("$x")
# If section starts w on A_DNA in state 1, calculate duration
    start_time = time_points[1];
    end_time = time_points[end];
    if states[1] == 1 && states[end] ==0
        first_transition = transitions_1_to_0[1] - start_time
        push!(on_state_durations, first_transition)
    # Loop for when the segment starts in the 1 state: Iterate over each transition from 1 to 0
        for i in 1:length(transitions_0_to_1)
        # Check if the submatrix starts in the 1 state and it's the first transition from 0 to 1
            if states[1] == 1 && i == 1
        # Subtract the first time point of 1 to 0 transition from the first time point of 0 to 1 transition
            duration = (transitions_1_to_0[i+1] - transitions_0_to_1[i])
            else
        # Otherwise, subtract the next consecutive time point of 1 to 0 transition from the current time point of 0 to 1 transition
            duration = (transitions_1_to_0[i+1] - transitions_0_to_1[i])
            end
        # Store the duration
            push!(on_state_durations, duration)
        end

    elseif states[end] == 1 && states[1] == 0 #segments starts in 0 state but ends in 1 state
        for i in 1:length(transitions_1_to_0)
            duration = transitions_1_to_0[i] - transitions_0_to_1[i]
            push!(on_state_durations, duration)
        end
        last_transition = end_time - transitions_0_to_1[end]
        push!(on_state_durations, last_transition)

    elseif states[1] == 1 && states[end] == 1 && length(transition_indices_0_to_1) > 0 #segment starts in 1 state and ends in 1 state
        first_transition = transitions_1_to_0[1] - start_time
        push!(on_state_durations, first_transition)
        for i in 1:length(transitions_0_to_1)-1
        # Check if the submatrix starts in the 1 state and it's the first transition from 0 to 1
            if states[1] == 1 && i == 1
        # Subtract the first time point of 1 to 0 transition from the first time point of 0 to 1 transition
            duration = (transitions_1_to_0[i+1] - transitions_0_to_1[i])
            else
        # Otherwise, subtract the next consecutive time point of 1 to 0 transition from the current time point of 0 to 1 transition
            duration = (transitions_1_to_0[i+1] - transitions_0_to_1[i])
            end
        # Store the duration
            push!(on_state_durations, duration)
        end
        last_transition = end_time - transitions_0_to_1[end]
        push!(on_state_durations, last_transition)
    elseif states[1]==0 && states[end]==0 && length(transition_indices_0_to_1) > 0 #Segment starts in 0 state and ends in 0 state
        for i in 1:length(transitions_0_to_1)
            duration = transitions_1_to_0[i] - transitions_0_to_1[i]
            push!(on_state_durations, duration)
        end
    elseif length(transition_indices_0_to_1)==0 && length(transition_indices_1_to_0)==0 && states[1]==1 && states[end]==1
        duration = end_time-start_time
        push!(on_state_durations, duration)
    elseif length(transition_indices_0_to_1)==0 && length(transition_indices_1_to_0)==0 && states[1]==0 && states[end]==0
        duration = 0
        push!(on_state_durations, duration)
    end
    total_on_state_duration = sum(on_state_durations)
    push!(total_on_state_durations, total_on_state_duration)
    
end
total_on_state_duration = round.(total_on_state_durations; digits = 6)
#On time per section
Ton_vs_Ttot = round.((total_on_state_durations / bin_width);digits = 6)
print(Ton_vs_Ttot)
plot(Ton_vs_Ttot)
# Calculate mean and mode
mean_Ton_vs_Ttot = round(mean(Ton_vs_Ttot),digits=6)
mode_Ton_vs_Ttot = round(StatsBase.mode(Ton_vs_Ttot),digits=6)

# Plot histogram
histogram(Ton_vs_Ttot, xlabel="Ton / Ttot", ylabel="Frequency", bins=:20, title="Histogram of Ton / Ttot",
    label=false);
vline!([mean_Ton_vs_Ttot], label="Mean", color=:red, linewidth=2, linestyle=:dash);
vline!([mode_Ton_vs_Ttot], label="Mode", color=:green, linewidth=2, linestyle=:dash);

# Print mean and mode values on the histogram
annotate!([(0.9, 40, text("Mean: $mean_Ton_vs_Ttot", :red)),
            (0.9, 37, text("Mode: $mode_Ton_vs_Ttot", :green))]);

plot!()

#analyze total durations on
mean_duration = []
mode_duration = []
std_duration = []
# Calculate mean
mean_duration = mean(total_on_state_durations);
# Calculate mode (using StatsBase package)
mode_duration = StatsBase.mode(total_on_state_durations);
# Calculate standard deviation
std_duration = std(total_on_state_durations);
# Define the number of significant figures you want to display
sig_figs = 5;
# Round the values to the specified number of significant figures
mean_duration_rounded = round(mean_duration; digits=sig_figs);
mode_duration_rounded = round(mode_duration; digits=sig_figs);
std_duration_rounded = round(std_duration; digits=sig_figs);
histogram(total_on_state_durations, bins=50, xlabel="Time promoter is in On state [Hrs]", ylabel="Frequency",
 title="Histogram of On State Durations (5 minute intervals)\nFeedback Removed, Induced Steady State", label= false, size = (800,600), xticks=10,
 xlims=(0.0,0.3333));


# Overlay summary statistics
vline!([mean_duration], label="Mean", linewidth=3, color=:red, linestyle=:dash);
vline!([mode_duration], label="Mode", linewidth=3, color=:purple, linestyle=:dash);


# Print values on the plot
annotate!([(0.26, 85, text("Mean: $mean_duration_rounded", :red, :left)),
            (0.26, 75, text("Mode: $mode_duration_rounded", :purple, :left)),
            (0.255, 65, text("Std Dev: $std_duration_rounded", :blue, :left))]);  

# Add legend
plot!() # This is needed to display the legend# This is needed to display the legend

#analyze total durations on
df_sums = DataFrame(sums = total_on_state_durations)
print(df_sums)
plot(total_on_state_durations, legend=false)
print(minimum(total_on_state_durations))
print(maximum(total_on_state_durations))

