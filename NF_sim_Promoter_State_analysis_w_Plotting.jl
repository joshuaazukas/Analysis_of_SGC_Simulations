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
p = [kOn => 10, kOff => 10, kOnt=> 0.000001, kOfft=> 100000, k => 3, kT => 1.5, deg_R => 0.03, deg_G => 0.008];
@named rn = ReactionSystem(rxs, t, [A, DNA, A_DNA, DNA_T, A_DNA_T, RNA, GFP], [kOn, kOff, kOnt, kOfft, k, kT, deg_R, deg_G]);
kOn = 10;
kOff = 10;
kOnt= 0.000001;
kOfft= 100000;
k = 3;
kT = 1.5;
deg_R = 0.03;
deg_G = 0.008;
num_simulations = 3;
T_of_sims = [];
A_DNA_sims = [];
GFP_sims = [];
T_A_DNA = [];
on_state_durations = [];
on_state_durations_sims = [];
total_on_state_duration = [];
total_on_state_durations = [];
total_on_state_durations_sims = [];
sub_matrices = [];
sub_matrix =[];
T_A_DNA =[];
bin_width = 0.3333; #also change in line 73
bin_width_min = Int(round.(bin_width * 60, digits = 1))
@time begin
    @threads for i in 1:num_simulations;

        T_A_DNA =[];
        dprob = DiscreteProblem(rn, u0, tspan, p)
        jprob = JumpProblem(rn, dprob, Direct())

        sol = solve(jprob, SSAStepper())
        
        T_A_DNA = (sol.t, sol[3,:]);

        T_GPF = (Float32(sol.t), Float32(sol[7,:]));

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
        bin_width = 0.3333;
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
        push!(GFP_sims, T_GPF)
        print("$i\n")
    end
end
Ton_vs_Ttot_per_simulation = []
length(total_on_state_durations_sims)

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
df_combined_Ton_vs_Ttot = DataFrame(comb = combined_Ton_vs_Ttot);
CSV.write("c://Users//jrazu//Desktop//JL//Promoter Statistics//combined_Ton_vs_Ttot$bin_width_min $kOn.csv", df_combined_Ton_vs_Ttot);
# Calculate mean and mode
mean_Ton_vs_Ttot_sims = round(mean(combined_Ton_vs_Ttot),digits=6);
mode_Ton_vs_Ttot_sims = round(StatsBase.mode(combined_Ton_vs_Ttot),digits=6);
stddev_Ton_vs_Ttot_sims = round(std(combined_Ton_vs_Ttot), digits=6);



# Plot histogram
Plots.histogram(combined_Ton_vs_Ttot, xlabel="Ton / Ttot", ylabel="PDF", bins=:20, title="Histogram of Ton / Ttot\n $bin_width_min min intervals ($num_simulations sims)",
    label=false, size = (800,800), legend=false, normalize=:probability);
vline!([mean_Ton_vs_Ttot_sims], label="Mean", color=:red, linewidth=2, linestyle=:dash);
vline!([mode_Ton_vs_Ttot_sims], label="Mode", color=:purple, linewidth=2, linestyle=:dash);

# Print mean and mode values on the histogram
annotate!([(0.02, 0.35, text("Mean: $mean_Ton_vs_Ttot_sims", :red)),
            (0.02, 0.33, text("Mode: $mode_Ton_vs_Ttot_sims", :purple)),
            (0.02, 0.31, text("Std: $stddev_Ton_vs_Ttot_sims", :blue)),
            (0.02, 0.28, text("kOn = $kOn, kOff = $kOff\nkOfft = $kOfft\n kOnt = $kOnt"))]);

plot!()
savefig("c://Users//jrazu//Desktop//JL//Promoter Statistics//combined_Ton_vs_Ttot$bin_width_min $kOn.png")
