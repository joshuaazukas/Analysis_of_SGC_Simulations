#evaluate parameter effect on promotor on state statistics
using Catalyst, DifferentialEquations, Plots, Interpolations, Statistics, Distributions, StatsBase, CSV, DataFrames

@parameters kOn kOff kOnt kOfft k kT deg_R deg_G
@variables t
@species A(t) DNA(t) A_DNA(t) DNA_T(t) A_DNA_T(t) RNA(t) GFP(t)

# Define the reaction network
rxs = [
    (@reaction kOn, A + DNA --> A_DNA),
    (@reaction kOff, A_DNA --> A + DNA),
    (@reaction kOnt, DNA + GFP --> DNA_T),
    (@reaction kOfft, DNA_T --> DNA + GFP),
    (@reaction kOnt, A_DNA + GFP --> A_DNA_T),
    (@reaction kOfft, A_DNA_T --> A_DNA),
    (@reaction k, A_DNA --> A_DNA + RNA),
    (@reaction kT, RNA --> RNA + GFP),
    (@reaction deg_R, RNA --> 0),
    (@reaction deg_G, GFP --> 0)
]

tspan = (0.0, 72)  # reaction time span
u0 = [A => 1, DNA => 1, A_DNA => 0, DNA_T => 0, A_DNA_T => 0, RNA => 1000, GFP => 950000]  # starting conditions
p_base = [1.0, 1.0, 0.0, 0.0, 3.0, 1.5, 0.03, 0.008]

@named rn = ReactionSystem(rxs, t, [A, DNA, A_DNA, DNA_T, A_DNA_T, RNA, GFP], [kOn, kOff, kOnt, kOfft, k, kT, deg_R, deg_G])

function run_simulation(rn, u0, tspan, p, num_runs)
    fractions_on_state = Float32[]
    total_changes_all = Float32[]
    frequency_changes_all = Float32[]

    for _ in 1:num_runs
        dprob = DiscreteProblem(rn, u0, tspan, p)
        jprob = JumpProblem(rn, dprob, Direct())
        sol = solve(jprob, SSAStepper())

        A_DNA = sol[3,:]

        # Analyze state duration
        on_durations_minutes = Float32[]
        in_1_state = false
        start_time = NaN

        for (i, val) in enumerate(A_DNA)
            if val == 1 && !in_1_state
                in_1_state = true
                start_time = sol.t[i]
            elseif val == 0 && in_1_state
                duration_minutes = (sol.t[i] - start_time) * 60
                push!(on_durations_minutes, duration_minutes)
                in_1_state = false
            end
        end

        if in_1_state
            duration_minutes = (tspan[2] - start_time) * 60
            push!(on_durations_minutes, duration_minutes)
        end

        total_on_time_minutes = sum(on_durations_minutes)
        total_simulation_time_minutes = (tspan[2] - tspan[1]) * 60
        fraction_on_time = total_on_time_minutes / total_simulation_time_minutes
        push!(fractions_on_state, fraction_on_time)

        # Analyze state changes
        transition_0_to_1 = 0
        transition_1_to_0 = 0

        for i in 2:length(A_DNA)
            if A_DNA[i-1] == 0 && A_DNA[i] == 1
                transition_0_to_1 += 1
            elseif A_DNA[i-1] == 1 && A_DNA[i] == 0
                transition_1_to_0 += 1
            end
        end

        total_changes = transition_0_to_1 + transition_1_to_0
        push!(total_changes_all, total_changes)

        total_time_hours = tspan[2] - tspan[1]
        average_frequency_per_hour = total_changes / total_time_hours
        push!(frequency_changes_all, average_frequency_per_hour)
    end

    avg_fraction_on_state = mean(fractions_on_state)
    std_fraction_on_state = std(fractions_on_state)
    avg_total_changes = mean(total_changes_all)
    std_total_changes = std(total_changes_all)
    avg_frequency_changes_per_hour = mean(frequency_changes_all)
    std_frequency_changes_per_hour = std(frequency_changes_all)

    return (avg_fraction_on_state, std_fraction_on_state, avg_total_changes, std_total_changes, avg_frequency_changes_per_hour, std_frequency_changes_per_hour)
end

num_runs = 2
multipliers = [10, 100, 1000, 10000, 10e5,10e6]

# Mapping parameter names to their indices
param_indices = Dict(:kOn => 1, :kOff => 2, :kOnt => 3, :kOfft => 4, :k => 5, :kT => 6, :deg_R => 7, :deg_G => 8)

# Store results for plotting
results = []

for multiplier in multipliers
    p = copy(p_base)
    p[param_indices[:kOn]] *= multiplier
    p[param_indices[:kOff]] *= multiplier
    (avg_fraction_on_state, std_fraction_on_state, avg_total_changes, std_total_changes, avg_frequency_changes_per_hour, std_frequency_changes_per_hour) = run_simulation(rn, u0, tspan, p, num_runs)

    println("Parameters with kOn and kOff multiplied by $multiplier:")
    println("Average fraction of time promoter is in On state: ", avg_fraction_on_state)
    println("Average total number of state changes: ", avg_total_changes)
    println("Average frequency of changes per hour: ", avg_frequency_changes_per_hour)

    push!(results, (multiplier, avg_fraction_on_state, std_fraction_on_state, avg_total_changes, std_total_changes, avg_frequency_changes_per_hour, std_frequency_changes_per_hour))
end

using Printf

# Plotting
multiplier_values = [result[1] for result in results]
avg_fractions = [result[2] for result in results]
std_fractions = [result[3] for result in results]
avg_frequencies = [result[6] for result in results]
std_frequencies = [result[7] for result in results]

# Custom tick labels for x-axis
x_tick_labels = [Printf.@sprintf("%.1f", log10(multiplier)) for multiplier in multiplier_values]

bar1 = bar(log.(multiplier_values), (avg_fractions), yerr=std_fractions, xlabel="kOn = kOff (log10)", ylabel="Average Fraction On State", title="Average Fraction On State vs kOn = kOff", legend=false, xticks=(log.(multiplier_values), x_tick_labels))
display(bar1)

bar2 = bar(log.(multiplier_values), log.(avg_frequencies), yerr=std_frequencies, xlabel="kOn = kOff (log10)", ylabel="Average Frequency of Changes per Hour (log10)", title="Average Frequency of Changes\nper Hour vs kOn = kOff", legend=false, xticks=(log.(multiplier_values), x_tick_labels))
display(bar2)