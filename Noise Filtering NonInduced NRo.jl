using Catalyst, Plots, Interpolations, JumpProcesses
using Statistics, Distributions, StatsBase
using CSV, DataFrames, FFTW, DSP
using SymbolicIndexingInterface: parameter_values, state_values
using SciMLStructures: Tunable, replace, replace!
using JLD2, UnPack
using XLSX
using GLM
#Z:/Single Cell Imaging/Nonind NR Comp/NR_ssNon_Comp.xlsx
NInd_data = XLSX.readxlsx("Z:/Single Cell Imaging/Nonind NR Comp/NR_ssNon_Comp.xlsx");
NInd_GFP = NInd_data["#GFP"];
NInd_GFP = NInd_GFP[:];
NInd_GFP = NInd_GFP[2:end,1:end]/1000000;
time = round.(collect(0:215).*0.33334, digits=2);

plot(NInd_GFP,legend=:false)

freq=0.00083333
cols = size(NInd_GFP,2);
fft_result1 = fft(NInd_GFP[:,1]*1000000)
freqs1 = fftfreq(length(time), 0.00083333)
plot(freqs1[1:div(end,2)], abs.(fft_result1[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum")

lp = Lowpass(0.000085, fs=0.00083333)
dmeth = Butterworth(2)
filt_GFP1 = filtfilt(digitalfilter(lp,dmeth), NInd_GFP[:,1]*1000000)
plot(time,filt_GFP1, label="Filtered GFP",color=:black)
plot!(time,NInd_GFP[:,1]*1000000, label="Raw GFP", xlabel="Time (Hrs)", ylabel="# of GFP molecules", color=:green)
plot();

for i in 1:cols
    fft_result = fft(NInd_GFP[:,i])
    freqs = fftfreq(length(time), 0.00083333)
    plot!(freqs[1:div(end,2)], abs.(fft_result[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum")
end
vline!([0.000085],color=:black,legend=:false)
plot!(yscale=:log10,legend=:false)

#define lowpass filter cutoff and frequency
lp = Lowpass(0.000085, fs=0.00083333)
dmeth = Butterworth(2)
#apply filter to GFP time trace example to remove high frequency noise
filt_NInd_GFP_d = []
NInd_GFP_means = []
for i in 1:cols
    k1 = mean(NInd_GFP[:,i])
    push!(NInd_GFP_means,k1)
    filt_GFP = filtfilt(digitalfilter(lp,dmeth), NInd_GFP[:,i])
    push!(filt_NInd_GFP_d,filt_GFP)
end
plot(time,filt_NInd_GFP_d,legend=:false)

NInd_mCh = NInd_data["#mCh"];
NInd_mCh = NInd_mCh[:];
NInd_mCh = NInd_mCh[2:end,1:end]/1000000;
plot(time,NInd_mCh,legend=:false,xlabel="#mCh (10^6)")
plot();

lp = Lowpass(0.000085, fs=0.00083333) #0.000085
dmeth = Butterworth(2)
filt_mCh1 = filtfilt(digitalfilter(lp,dmeth), NInd_mCh[:,1]*1000000)
plot(time,filt_mCh1, label="Filtered mCh",color=:black)
plot!(time,NInd_mCh[:,1]*1000000, label="Raw mCh", xlabel="Time (Hrs)", ylabel="# of mCh molecules", color=:red)
plot();

for i in 1:cols
    fft_result = fft(NInd_mCh[:,i])
    freqs = fftfreq(length(time), 0.00083333)
    plot!(freqs[1:div(end,2)], abs.(fft_result[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum")
end
plot!(legend=:false, xlimit=(-0.00001,0.0001))

filt_NInd_mCh_d = []
NInd_mCh_means = []
for i in 1:cols
    k1 = mean(NInd_mCh[:,i])
    push!(NInd_mCh_means,k1)
    filt_mCh = filtfilt(digitalfilter(lp,dmeth), NInd_mCh[:,i])
    push!(filt_NInd_mCh_d,filt_mCh)
end
plot(time,filt_NInd_mCh_d,legend=:false)
filt_mCh_traces = Float64.(DataFrame([Symbol("Col$i") => vec for (i, vec) in enumerate(filt_NInd_mCh_d*1000000)]))
filt_GFP_traces = Float64.(DataFrame([Symbol("Col$i") => vec for (i, vec) in enumerate(filt_NInd_GFP_d*1000000)]))

CSV.write("Z:/Single Cell Imaging/Nonind NR Comp/filt_#mCh.csv", filt_mCh_traces)
CSV.write("Z:/Single Cell Imaging/Nonind NR Comp/filt_#GFP.csv", filt_GFP_traces)
