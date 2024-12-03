using Catalyst, Plots, Interpolations, JumpProcesses
using Statistics, Distributions, StatsBase
using CSV, DataFrames, FFTW, DSP
using SymbolicIndexingInterface: parameter_values, state_values
using SciMLStructures: Tunable, replace, replace!
using JLD2, UnPack
using XLSX
using GLM

BFCexample_GFP = XLSX.readxlsx("C:/Users/Strey Lab/Documents/GitHub/Analysis_of_SGC_Simulations/StreyCats/Broken Feedback Circuit/Broken_Feedback_Circuit_High_Expression_Example.xlsx");
GFPexample = BFCexample_GFP["Sheet1"];
GFPexample = GFPexample[:];
example_GFP = GFPexample[2:end,2:end]/1000000
plot(example_GFP)
time = round.(GFPexample[2:end,1].*0.33334, digits=2)
freq=0.00083333
lp = Lowpass(0.00006, fs=0.00083333)
dmeth = Butterworth(2)
filt_GFP_example = filtfilt(digitalfilter(lp,dmeth), example_GFP)
plot(time,filt_GFP_example)
plot!(time,example_GFP)
Example_Filt_minus=example_GFP-filt_GFP_example
plot(Example_Filt_minus)
std_Exam_Filt_min = std(Example_Filt_minus)
ndist = std_Exam_Filt_min.*randn(216)
histogram(ndist, bins=10)
filt_w_ndist = filt_GFP_example + ndist
plot(time,filt_w_ndist)
plot!(time,filt_GFP_example)
plot!(time,example_GFP)

plot(Example_Filt_minus)

# Perform FFT
fft_result = fft(Example_Filt_minus)
freqs = fftfreq(length(time), 1/0.00083333)

# Plot the frequency spectrum
plot(freqs[1:div(end,2)], abs.(fft_result[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum")


BFC_GFP = XLSX.readxlsx("C:/Users/Strey Lab/Documents/GitHub/Analysis_of_SGC_Simulations/StreyCats/Broken Feedback Circuit/Broken Feedback Circuit Time Traces 5_21_24 GFP.xlsx");
GFPex = BFC_GFP["Sheet1"];
GFPex = GFPex[:];
time = round.(GFPex[2:end,1].*0.33334, digits=2)
ex_GFP = GFPex[2:end,2:end]/1000000
plot(time,ex_GFP,legend=:false)
#GFP = GFP.*0.18
cols = size(ex_GFP,2);

freq=0.00083333
lp = Lowpass(0.00006, fs=0.00083333)
dmeth = Butterworth(2)
filt_GFP = filtfilt(digitalfilter(lp,dmeth), ex_GFP[:,1])

filt_GFP_traces = []
filt_GFP_noise_traces = []
GFP_means = []
GFP_STD = []

for i in 1:cols
    filt_GFP=[]
    stdGFP=[]
    filt_noise=[]
    k1 = mean(ex_GFP[:,i])
    push!(GFP_means,k1)
    filt_GFP = filtfilt(digitalfilter(lp,dmeth), ex_GFP[:,i])
    push!(filt_GFP_traces,filt_GFP)
    Ex_Filt=ex_GFP[:,i]-filt_GFP
    stdGFP=std(Ex_Filt)
    push!(GFP_STD,stdGFP)
    noise=stdGFP.*randn(216)
    filt_noise = filt_GFP+noise
    push!(filt_GFP_noise_traces, filt_noise)
end

GFP_exp_temp_u = reduce(vcat, transpose.(GFP_means))
GFP_noise_std = reduce(vcat, transpose.(GFP_STD))

plot(filt_GFP_noise_traces, legend=:false)
scatter(GFP_exp_temp_u,GFP_noise_std,xlabel="Temporal Mean Experimental Traces (# GFP x10^6)", ylabel="standard deviation of lowpass filtered\ndata subtracted experimental traces",legend=:false)

scatter(log10.(GFP_exp_temp_u),log10.(GFP_noise_std),xlabel="Temporal Mean Experimental Traces (# GFP x10^6)", ylabel="standard deviation of lowpass filtered\ndata subtracted experimental traces",legend=:false)

df = DataFrame(u=GFP_exp_temp_u,std=GFP_noise_std)

CSV.write("Temporal_mean_std_filt_data_min_exp.csv", df)