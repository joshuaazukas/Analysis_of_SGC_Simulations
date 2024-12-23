using Catalyst, Plots, Interpolations, JumpProcesses
using Statistics, Distributions, StatsBase
using CSV, DataFrames, FFTW, DSP
using SymbolicIndexingInterface: parameter_values, state_values
using SciMLStructures: Tunable, replace, replace!
using JLD2, UnPack
using XLSX
using GLM

#Get Experimental time trace data High Expression Example
HInd_data = XLSX.readxlsx("Z:/Single Cell Imaging/4_19_24_NRo_ss//NRo 500 Dox Induced Steady State 4_19_24 C1-3.xlsx");
HInd_GFP = HInd_data["#GFP"];
HInd_GFP = HInd_GFP[:];
HInd_GFP = HInd_GFP[2:end,1:end]/1000000;
time = round.(collect(0:215).*0.33334, digits=2);
plot(HInd_GFP,legend=:false)
HInd_GFP
#frequency (1/1200s) measurement every 20 minutes
freq=0.00083333
cols = size(HInd_GFP,2);
fft_result1 = fft(HInd_GFP[:,1]*1000000)
freqs1 = fftfreq(length(time), 0.00083333)
plot(freqs1[1:div(end,2)], abs.(fft_result1[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum")
lp = Lowpass(0.000085, fs=0.00083333)
dmeth = Butterworth(2)
filt_GFP1 = filtfilt(digitalfilter(lp,dmeth), HInd_GFP[:,1]*1000000)
plot(time,filt_GFP1, label="Filtered GFP",color=:black)
plot!(time,HInd_GFP[:,1]*1000000, label="Raw GFP", xlabel="Time (Hrs)", ylabel="# of GFP molecules", color=:green)
plot();
for i in 1:cols
    fft_result = fft(HInd_GFP[:,i])
    freqs = fftfreq(length(time), 0.00083333)
    plot!(freqs[1:div(end,2)], abs.(fft_result[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum")
end
vline!([0.000085],color=:black,legend=:false)
plot!(yscale=:log10,legend=:false)

#define lowpass filter cutoff and frequency
lp = Lowpass(0.000085, fs=0.00083333)
dmeth = Butterworth(2)
#apply filter to GFP time trace example to remove high frequency noise

filt_HInd_GFP_d = []
HInd_GFP_means = []
for i in 1:cols
    k1 = mean(HInd_GFP[:,i])
    push!(HInd_GFP_means,k1)
    filt_GFP = filtfilt(digitalfilter(lp,dmeth), HInd_GFP[:,i])
    push!(filt_HInd_GFP_d,filt_GFP)
end
plot(time,filt_HInd_GFP_d,legend=:false)

HInd_mCh = HInd_data["#mCh"];
HInd_mCh = HInd_mCh[:];
HInd_mCh = HInd_mCh[2:end,1:end]/1000000;
plot(time,HInd_mCh,legend=:false,xlabel="#mCh (10^6)")
plot();

lp = Lowpass(0.000001, fs=0.00083333) #0.000085
dmeth = Butterworth(2)
filt_mCh1 = filtfilt(digitalfilter(lp,dmeth), HInd_mCh[:,1]*1000000)
plot(time,filt_mCh1, label="Filtered mCh",color=:black)
plot!(time,HInd_mCh[:,1]*1000000, label="Raw mCh", xlabel="Time (Hrs)", ylabel="# of mCh molecules", color=:red)
plot()
for i in 1:cols
    fft_result = fft(HInd_mCh[:,i])
    freqs = fftfreq(length(time), 0.00083333)
    plot!(freqs[1:div(end,2)], abs.(fft_result[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum")
end
plot!(legend=:false, xlimit=(-0.00001,0.0001))

filt_HInd_mCh_d = []
HInd_mCh_means = []
for i in 1:cols
    k1 = mean(HInd_mCh[:,i])
    push!(HInd_mCh_means,k1)
    filt_mCh = filtfilt(digitalfilter(lp,dmeth), HInd_mCh[:,i])
    push!(filt_HInd_mCh_d,filt_mCh)
end
plot(time,filt_HInd_mCh_d,legend=:false)
filt_mCh_traces = Float64.(DataFrame([Symbol("Col$i") => vec for (i, vec) in enumerate(filt_HInd_mCh_d*1000000)]))
filt_GFP_traces = Float64.(DataFrame([Symbol("Col$i") => vec for (i, vec) in enumerate(filt_HInd_GFP_d*1000000)]))

CSV.write("Z:/Single Cell Imaging/4_19_24_NRo_ss//filt_#mCh.csv", filt_mCh_traces)
CSV.write("Z:/Single Cell Imaging/4_19_24_NRo_ss//filt_#GFP.csv", filt_GFP_traces)

plot(time,filt_mCh_traces[:,20]*1000000, label="mCh Raw",color=:red,linestyle=:dash)
plot!(time,HInd_mCh[:,20]*1000000, label="mCh Filt",color=:red)
plot!(time,HInd_GFP[:,20]*1000000, label="GFP Raw",color=:green)
plot!(time,filt_GFP_traces[:,20]*1000000, label="GFP filt",color=:green,linestyle=:dash, ylabel="# of molecules",xlabel="time (Hrs)")

plot(log10.(filt_mCh_traces[:,20]),log10.(filt_GFP_traces[:,20]),label="Filtered",xlabel="Log10 # of GFP molecules",ylabel="Log10 # of mCh molecules",color=:black)
plot!(log10.(HInd_mCh[:,20]*1000000),log10.(HInd_GFP[:,20]*1000000),label="Raw",color=:orange)

plot(filt_mCh_traces[:,1],filt_GFP_traces[:,1],label="Filtered",xlabel="# of GFP molecules",ylabel="# of mCh molecules",yscale=:log10,xscale=:log10)
plot!((HInd_mCh[:,1]*1000000),(HInd_GFP[:,1]*1000000),label="Raw",yscale=:log10,xscale=:log10)