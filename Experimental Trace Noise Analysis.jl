using Catalyst, Plots, Interpolations, JumpProcesses
using Statistics, Distributions, StatsBase
using CSV, DataFrames, FFTW, DSP
using SymbolicIndexingInterface: parameter_values, state_values
using SciMLStructures: Tunable, replace, replace!
using JLD2, UnPack
using XLSX
using GLM

#Get Experimental time trace data High Expression Example
BFCexample_GFP = XLSX.readxlsx("C:/Users/Strey Lab/Documents/GitHub/Analysis_of_SGC_Simulations/StreyCats/Broken Feedback Circuit/Broken_Feedback_Circuit_High_Expression_Example.xlsx");
GFPexample = BFCexample_GFP["Sheet1"];
GFPexample = GFPexample[:];
example_GFP = GFPexample[2:end,2:end]/1000000
plot(example_GFP)
time = round.(GFPexample[2:end,1].*0.33334, digits=2)
time
#frequency (1/1200s) measurement every 20 minutes
freq=0.00083333
#define lowpass filter cutoff and frequency
lp = Lowpass(0.00006, fs=0.00083333)
dmeth = Butterworth(2)
#apply filter to GFP time trace example to remove high frequency noise
filt_GFP_example = filtfilt(digitalfilter(lp,dmeth), example_GFP)
#Compare Original and filtered time traces
plot(time,filt_GFP_example)
plot!(time,example_GFP, label="exp")
#Subtract filtered trace from experimental trace to get high frequency noise
Example_Filt_minus=example_GFP-filt_GFP_example
plot(Example_Filt_minus)
#calculate the standard deviation of the high frequency noise and use it to create a random array of values around 0 with the same standard deviation
std_Exam_Filt_min = std(Example_Filt_minus)
ndist = std_Exam_Filt_min.*randn(216)
histogram(ndist, bins=10)
#add the random noise array to the filtered experimental time trace
filt_w_ndist = filt_GFP_example + ndist
#Compare the original experimental, filtered, and filtered with added noise traces
plot(time,filt_w_ndist, label="lowpass Filt + Noise")
plot!(time,filt_GFP_example, label="Lowpass Filt.")
plot!(time,example_GFP, label="exp")

# visualize the frequency spectrum of the noise
fft_result = fft(example_GFP)
freqs = fftfreq(length(time), 0.00083333)
plot(freqs[1:div(end,2)], abs.(fft_result[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum")

Bkg = XLSX.readxlsx("C:/Users/Strey Lab/Documents/GitHub/Analysis_of_SGC_Simulations/StreyCats/Background/Background Traces.xlsx");
Cors = Bkg["Con3_Cor"];
OGs = Bkg["OG #GFP"];
time = round.(collect(0:215).*0.33334, digits=2);
Cor = Cors[:];
OG = OGs[:];
#plot Raw Data
Cols = size(Cor,2);
Cols



#Load all experimental data and produce a plot
BFC_GFP = XLSX.readxlsx("C:/Users/Strey Lab/Documents/GitHub/Analysis_of_SGC_Simulations/StreyCats/Broken Feedback Circuit/Broken Feedback Circuit Time Traces 5_21_24 GFP.xlsx");
GFPex = BFC_GFP["Sheet1"];
GFPex = GFPex[:];
time = round.(GFPex[2:end,1].*0.33334, digits=2)
ex_GFP = GFPex[2:end,2:end]/1000000
plot(time,ex_GFP,legend=:false)
cols = size(ex_GFP,2);
ex_GFP_m = mean(ex_GFP, dims=2)
#calculate the temporal standard deviation for prefiltered experimental time traces
exp_temp_std = []
for i in 1:cols
    exp_stdGFP=std(ex_GFP[:,i])
    push!(exp_temp_std,exp_stdGFP)
end
histogram(exp_temp_std, bins=30, normed=:probability)
plot();
fft_res_list = [];
for i in 1:cols
    fft_result = fft(ex_GFP[:,i])
    freqs = fftfreq(length(time), 0.00083333)
    push!(fft_res_list, fft_result)
    plot!(freqs[1:div(end,2)], abs.(fft_result[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum")
end
plot!(yscale=:log10,legend=:false)

fft_res_l = reduce(hcat,fft_res_list)
fft_res_m = mean(fft_res_l,dims=2)
plot(freqs[1:div(end,2)], abs.(fft_res_m[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum",yscale=:log10)
#determine the number of columns/samples in the data
cols = size(ex_GFP,2);
#define frequency and filter (repeated from above)
freq=0.00083333
lp = Lowpass(0.00006, fs=0.00083333)
dmeth = Butterworth(2)
filt_GFP = filtfilt(digitalfilter(lp,dmeth), ex_GFP[:,1])

#set up arrays to store outputs from the loop below
filt_GFP_traces = []
filt_GFP_noise_traces = []
GFP_means = []
filtGFP_STD=[]
GFP_STD = []
# loop for performing the same analysis and noise addition as performed on the example time trace above
# stores: Temporal mean of each trace prior to filtering, Filtered experimental traces, Standard deviation of the noise, and the filtered traces with the added noise
for i in 1:cols
    filt_GFP=[]
    stdGFP=[]
    filt_noise=[]
    Ex_Filt=[]
    stdfilt_GFP = []
    k1 = mean(ex_GFP[:,i])
    push!(GFP_means,k1)
    filt_GFP = filtfilt(digitalfilter(lp,dmeth), ex_GFP[:,i])
    push!(filt_GFP_traces,filt_GFP)
    stdfilt_GFP = std(filt_GFP)
    push!(filtGFP_STD,stdfilt_GFP)
    Ex_Filt=ex_GFP[:,i]-filt_GFP
    stdGFP=std(Ex_Filt)
    push!(GFP_STD,stdGFP)
    noise=stdGFP.*randn(216)
    filt_noise = filt_GFP+noise
    push!(filt_GFP_noise_traces, filt_noise)
end
plot(filt_GFP_traces,legend=:false)
# store the temporal mean and standard deviation of the noise for each time trace
GFP_exp_temp_u = reduce(vcat, transpose.(GFP_means))
GFP_noise_std = reduce(vcat, transpose.(GFP_STD))
#Store the temporal standard deviation of the filtered experimental time traces (befor adding noise back)
filt_temp_std = reduce(vcat, transpose.(filtGFP_STD))
#visualize the filtered experimental data with noise added back
plot(filt_GFP_noise_traces, legend=:false)
#format the filtered + noise traces into data frame
dffilt_GFP_noise_traces = DataFrame([Symbol("Col$i") => vec for (i, vec) in enumerate(filt_GFP_noise_traces)])
#calculate the temporal standard deviation for filtered + noise time traces and compare the distribution to prefiltered data
exp_filt_noise_temp_std = []
for i in 1:cols
    exp_filt_noise_stdGFP=std(dffilt_GFP_noise_traces[:,i])
    push!(exp_filt_noise_temp_std,exp_filt_noise_stdGFP)
end
histogram(exp_filt_noise_temp_std, bins=30, normed=:probability,color="grey",alpha=0.5, label="Filtered + Noise Traces")
histogram!(exp_temp_std, bins=30, normed=:probability,color="red",alpha=0.5,label="Raw Experimental Traces",xlabel="Temporal Std",ylabel="Probability")
#compare distributions of temporal standard deviation of the experimental traces and filtered experimental traces
histogram(filt_temp_std, bins=30, normed=:probability,color="grey",alpha=0.5, label="Filtered Traces")
histogram!(exp_temp_std, bins=30, normed=:probability,color="red",alpha=0.5,label="Raw Experimental Traces",xlabel="Temporal Std",ylabel="Probability")
#plot the standard deviation of the noise against the temporal mean of the experimental data
scatter(GFP_exp_temp_u,GFP_noise_std,xlabel="Temporal Mean Experimental Traces (# GFP x10^6)", ylabel="standard deviation of lowpass filtered\ndata subtracted experimental traces",legend=:false)
#plot the standard deviation of the noise against the log10 transformed temporal mean of the experimental data
scatter(log10.(GFP_exp_temp_u),log10.(GFP_noise_std),xlabel="Temporal Mean Experimental Traces (# GFP x10^6)", ylabel="standard deviation of lowpass filtered\ndata subtracted experimental traces",legend=:false)


#save the data as a 
#df = DataFrame(u=GFP_exp_temp_u,std=GFP_noise_std)
#CSV.write("Temporal_mean_std_filt_data_min_exp.csv", df)