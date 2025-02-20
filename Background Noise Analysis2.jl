using Catalyst, Plots, Interpolations, JumpProcesses
using Statistics, Distributions, StatsBase
using CSV, DataFrames, FFTW, DSP
using SymbolicIndexingInterface: parameter_values, state_values
using SciMLStructures: Tunable, replace, replace!
using JLD2, UnPack
using XLSX
using GLM
using DataFrames


# Linear detrending
function detrend_linear_columns_add_means(data::Matrix)
    # Initialize a matrix to store detrended data
    detrended_data = similar(data)
    detrended_datam = similar(data)
    # Loop through each column of the matrix
    for j in 1:size(data, 2)
        signal = Float64.(data[:, j])  # Extract the signal from column j
        time = 1:length(signal)  # Time vector (row indices)
        m = mean(signal,dims=1)
        # Create a DataFrame for GLM
        df = DataFrame(time=time, signal=signal)

        # Fit a linear model: signal ~ time
        lm_model = lm(@formula(signal ~ time), df)

        # Predict the trend
        trend = predict(lm_model)

        # Subtract the trend to detrend the signal
        detrended_data[:, j] = signal .- trend
        detrended_datam[:,j] = m .+ detrended_data[:, j]
    end
    return detrended_datam
end

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

selected_indices = sample(1:size(ex_GFP, 2), 18, replace=false)
selected_indices = [217, 826, 21, 373, 189, 348, 718, 669, 454, 433, 416, 320, 681, 472, 336, 877, 134, 893]
# Create a new matrix with the selected columns
Sel_ex_GFP = ex_GFP[:, selected_indices]

plot();
fft_res_list = [];
for i in 1:18
    mc = mean(Sel_ex_GFP[:,i])
    ex_GFPm = Sel_ex_GFP[:,i]./mc
    fft_result = fft(ex_GFPm)
    freqs = fftfreq(length(time), 0.00083333)
    push!(fft_res_list, fft_result)
    plot!(freqs[1:div(end,2)], abs.(fft_result[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum")
end
plot!(yscale=:log10,legend=:false)

fft_res_l = reduce(hcat,fft_res_list)
fft_res_m = mean(fft_res_l,dims=2)
plot(freqs[1:div(end,2)], abs.(fft_res_m[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum",yscale=:log10)

Bkg = XLSX.readxlsx("C:/Users/Strey Lab/Documents/GitHub/Analysis_of_SGC_Simulations/StreyCats/Background/Background Traces.xlsx");
Cors = Bkg["Cor #GFP"];
OGs = Bkg["OG #GFP"];
time = round.(collect(0:215).*0.33334, digits=2);
Cor = Cors[:];
OG = OGs[:];
#plot Raw Data
Cols = size(Cor,2);
Cols
plot();
for i in 1:18
    plot!(time,Cor[:,i]./mean(Cor[:,i]))
end
plot!()

plot(time,OG[:,1])
plot(time, Cor[:,1]./mean(Cor[:,1]), label="Cor")
plot!(time, OG[:,1]./mean(OG[:,1]), label="OG")
detrended_datam = detrend_linear_columns_add_means(OG)

plot();
for i in 1:18
    plot!(time,detrended_datam[:,i]/1000000)
end
plot!(title=:"Linear Detrended Data + Mean", xlabel="Time (Hrs)", ylabel = "# GFP Molecules (10^6)",legend=:false)

plot();
bkg_fft_res = [];
for i in 1:18
    m = mean(OG[:,i])
    Cor1 = OG[:,i]./m
    fft_result = fft(Cor1)
    freqs = fftfreq(length(time),0.00083333)
    plot!(freqs[1:div(end,2)], abs.(fft_result[1:div(end,2)]), title="Fourier Tansform", xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum")
    push!(bkg_fft_res, fft_result)
    #plot!(freqs, abs.(fft_result), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum",color=:red,xlim=(0,0.00045))
end
plot!(legend=:false,yscale=:log10)
bkg_fft_res_l = reduce(hcat,bkg_fft_res)
bkg_fft_res_m = mean(bkg_fft_res_l,dims=2)
plot(freqs[1:div(end,2)], abs.(bkg_fft_res_m[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Mean Bkg Frequency Spectrum",color=:red,yscale=:log10)
plot!(freqs[1:div(end,2)], abs.(fft_res_m[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Mean Cell Frequency Spectrum",color=:blue,title="Fourier Transform of Cell and Background Fluor.",yscale=:log10)
vline!([0.000085], color=:black, linestyle=:dash, label="LP: 8.5x10^-5 Hz")
plot();

bkg_fft_res = [];
for i in 1:18
    m = mean(detrended_datam[:,i])
    detrended_datamm = detrended_datam[:,i]./m
    fft_result = fft(Float64.(detrended_datamm))
    magnitude = abs.(fft_result) ./ 216
    freqs = fftfreq(length(time),0.00083333)
    positive_indices = findall(>=(0), freqs)
    positive_freqs = freqs[positive_indices]
    positive_magnitude = magnitude[positive_indices]
    plot!(positive_freqs, positive_magnitude, title="Fourier Tansform", xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum")
    push!(bkg_fft_res, fft_result)
    #plot!(freqs, abs.(fft_result), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum",color=:red,xlim=(0,0.00045))
end
plot!(legend=:false,yscale=:log10)
bkg_fft_res_l = reduce(hcat,bkg_fft_res)

bkg_fft_res_m = mean(bkg_fft_res_l,dims=2)
plot(freqs[1:div(end,2)], abs.(bkg_fft_res_m[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Mean Bkg Frequency Spectrum",color=:red,yscale=:log10)
plot!(freqs[1:div(end,2)], abs.(fft_res_m[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Mean Cell Frequency Spectrum",color=:blue,title="Fourier Transform of Cell and Background Fluor.",yscale=:log10)
vline!([0.00009], color=:black, linestyle=:dash, label="LP: 8.5x10^-5 Hz")
savefig(joinpath("StreyCats/Background/FTCellBkg.png"))
plot();
fft_res_list = [];
for i in 1:cols
    mc = mean(ex_GFP[:,i])
    ex_GFPm = ex_GFP[:,i]./mc
    fft_result = fft(ex_GFPm)
    freqs = fftfreq(length(time), 0.00083333)
    push!(fft_res_list, fft_result)
    plot!(freqs[1:div(end,2)], abs.(fft_result[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum",color=:black)
end
plot!(yscale=:log10,legend=:false)
bkg_fft_res = [];
for i in 1:18
    m = mean(detrended_datam[:,i])
    detrended_datamm = detrended_datam[:,i]./m
    fft_result = fft(Float64.(detrended_datamm))
    magnitude = abs.(fft_result) ./ 216
    freqs = fftfreq(length(time),0.00083333)
    positive_indices = findall(>=(0), freqs)
    positive_freqs = freqs[positive_indices]
    positive_magnitude = magnitude[positive_indices]
    plot!(positive_freqs, positive_magnitude, title="Fourier Tansform", xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum",color=:red)
    push!(bkg_fft_res, fft_result)
    #plot!(freqs, abs.(fft_result), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum",color=:red,xlim=(0,0.00045))
end
plot!(legend=:false,yscale=:log10,color=:red)

#define frequency and filter (repeated from above)
freq=0.00083333
lp = Lowpass(0.00009, fs=0.00083333)
dmeth = Butterworth(2)
filt_GFP = filtfilt(digitalfilter(lp,dmeth), ex_GFP[:,1])

#set up arrays to store outputs from the loop below
filt_GFP1=[]
filt_GFP_traces1 = []
filt_GFP_noise_traces1 = []
GFP_means1 = []
filtGFP_STD1=[]
GFP_STD1 = []
# loop for performing the same analysis and noise addition as performed on the example time trace above
# stores: Temporal mean of each trace prior to filtering, Filtered experimental traces, Standard deviation of the noise, and the filtered traces with the added noise
for i in 1:cols
    filt_GFP1=[]
    stdGFP1=[]
    filt_noise1=[]
    Ex_Filt1=[]
    stdfilt_GFP1= []
    k1 = mean(ex_GFP[:,i])
    push!(GFP_means1,k1)
    filt_GFP = filtfilt(digitalfilter(lp,dmeth), ex_GFP[:,i])
    push!(filt_GFP_traces1,filt_GFP)
    stdfilt_GFP = std(filt_GFP)
    push!(filtGFP_STD1,stdfilt_GFP)
    Ex_Filt=ex_GFP[:,i]-filt_GFP
    stdGFP=std(Ex_Filt)
    push!(GFP_STD1,stdGFP)
    noise=stdGFP.*randn(216)
    filt_noise = filt_GFP+noise
    push!(filt_GFP_noise_traces1, filt_noise)
end
plot(filt_GFP_traces1,legend=:false)
scatter(GFP_means1,GFP_STD1)
plot(time, filt_GFP_traces1, label="Filtered ")

#define frequency and filter (repeated from above)
freq=0.00083333
lp = Lowpass(0.000085, fs=0.00083333)
dmeth = Butterworth(2)
filt_GFP = filtfilt(digitalfilter(lp,dmeth), ex_GFP[:,1])

#set up arrays to store outputs from the loop below
filt_GFP2=[]
filt_GFP_traces2 = []
filt_GFP_noise_traces2 = []
GFP_means2 = []
filtGFP_STD2=[]
GFP_STD2 = []
# loop for performing the same analysis and noise addition as performed on the example time trace above
# stores: Temporal mean of each trace prior to filtering, Filtered experimental traces, Standard deviation of the noise, and the filtered traces with the added noise
for i in 1:cols
    filt_GFP2=[]
    stdGFP2=[]
    filt_noise2=[]
    Ex_Filt2=[]
    stdfilt_GFP2= []
    k1 = mean(ex_GFP[:,i])
    push!(GFP_means2,k1)
    filt_GFP = filtfilt(digitalfilter(lp,dmeth), ex_GFP[:,i])
    push!(filt_GFP_traces2,filt_GFP)
    stdfilt_GFP = std(filt_GFP)
    push!(filtGFP_STD2,stdfilt_GFP)
    Ex_Filt=ex_GFP[:,i]-filt_GFP
    stdGFP=std(Ex_Filt)
    push!(GFP_STD2,stdGFP)
    noise=stdGFP.*randn(216)
    filt_noise = filt_GFP+noise
    push!(filt_GFP_noise_traces2, filt_noise)
end
plot(filt_GFP_traces2,legend=:false)
scatter(GFP_means1,GFP_STD1, label=:"LP:0.00009", color=:red)
scatter!(GFP_means2,GFP_STD2, label=:"LP:0.000085",color=:blue)


