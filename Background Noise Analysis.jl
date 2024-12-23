using Catalyst, Plots, Interpolations, JumpProcesses
using Statistics, Distributions, StatsBase
using CSV, DataFrames, FFTW, DSP
using SymbolicIndexingInterface: parameter_values, state_values
using SciMLStructures: Tunable, replace, replace!
using JLD2, UnPack
using XLSX
using GLM
using DataFrames
function movmean(data::Matrix{T}, window_size::Int) where T
    n = length(data)
    smoothed = similar(data)
    half_window = div(window_size, 2)
    for i in 1:n
        start_idx = max(1, i - half_window)
        end_idx = min(n, i + half_window)
        smoothed[i] = mean(data[start_idx:end_idx])
    end
    return smoothed
end

Bkg = XLSX.readxlsx("C:/Users/Strey Lab/Documents/GitHub/Analysis_of_SGC_Simulations/StreyCats/Background/Background Traces.xlsx");
Cors = Bkg["Con3_Cor"];
OGs = Bkg["OG #GFP"];
time = round.(collect(0:215).*0.33334, digits=2);
Cor = Cors[:];
OG = OGs[:];
#plot Raw Data
Cols = size(Cor,2);
Cols
plot();
for i in 1:Cols
    #plot!(Cor[:,i],legend=:false)
    plot!(time,OG[:,i]/1000000,legend=:false)
end
plot!(title="Raw Background Traces (n=18)", ylabel="# of GFP Molecules (x 10^6)", xlabel="Time (Hrs)")

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

# Linear detrending
function detrend_linear_columns(data::Matrix)
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

    end
    return detrended_data
    
end

detrended_data = detrend_linear_columns(OG)
detrended_datam = detrend_linear_columns_add_means(OG)

plot();
for i in 1:18
    plot!(time,detrended_data[:,i])
end
plot!(title=:"Linear Detrended Data",legend=:false)

plot();
for i in 1:18
    plot!(time,detrended_datam[:,i]/1000000)
end
plot!(title=:"Linear Detrended Data + Mean", xlabel="Time (Hrs)", ylabel = "# GFP Molecules (10^6)",legend=:false)

det_datam = Float64.(DataFrame(detrended_datam, :auto))
CSV.write("C:/Users/Strey Lab/Documents/GitHub/Analysis_of_SGC_Simulations/StreyCats/detrended_BKG_Traces.csv", det_datam)
ms=[];
stds=[];
for i in 1:18
    m = mean(detrended_datam[:,i])
    stdd = std(detrended_datam[:,i])
    push!(stds,stdd)
    push!(ms, m)
end
mstd = mean(stds)
mm = mean(ms)
df = DataFrame(u=ms,std=stds)
noise=mstd.*randn(216)
filt_noise = mm+noise # need to make mm a vector 216 points long

CSV.write("C:/Users/Strey Lab/Documents/GitHub/Analysis_of_SGC_Simulations/StreyCats/Background/detrended_BKG_Traces_u_std.csv", df)
scatter(ms,stds, xlabel="Mean of detrended + Mean Traces", ylabel="Standard Deviation")
#subtract mean from data traces

OGms = mean(OG, dims=1)
OGmsub = OG .-OGms
plot(OGmsub)
plot(time, OGmsub./([OGmsub[1]]))
# Moving Average detrending
window_size = 50 
smoothed = movmean(OGms, window_size)
plot(smoothed)
OGm_det = OGms - smoothed
plot(OGm_det)
OG_det_d = []
for i  in 1:18
    smoothed = movmean(OG[:,i],window_size)
    OG_det = OG[:,i] - smoothed
    push!(OG_det_d,OG_det)
end
OG_det_d_mat = Float64.(DataFrame([Symbol("Col$i") => vec for (i, vec) in enumerate(OG_det_d)]))
plot(OG_det_d, legend=:false)

sr = 0.00083333
per = [];
freqsp_d = [];
#FFT
plot();
for i in 1:18
    fft_result = fft(Float64.(detrended_datam[:,i]))
    magnitude = abs.(fft_result) ./ 216
    freqs = fftfreq(length(time),0.00083333)
    positive_indices = findall(>=(0), freqs)
    positive_freqs = freqs[positive_indices]
    positive_magnitude = magnitude[positive_indices]
    plot!(positive_freqs, positive_magnitude, title="Fourier Tansform", xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum")
    #plot!(freqs, abs.(fft_result), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum",color=:red,xlim=(0,0.00045))
end
plot!(legend=:false)
vline!([0.000085],color=:black,legend=:false)

fft_result=[];
magnitude=[];
positive_magnitude=[];
positive_freqs=[];
plot();
for i in 1:18
    fft_result = fft(Float64.(detrended_data[:,i]))
    magnitude = abs.(fft_result) ./ 216
    
    freqs = fftfreq(length(time),0.00083333)
    positive_indices = findall(>=(0), freqs)
    positive_freqs = freqs[positive_indices]
    positive_magnitude = magnitude[positive_indices]
    
    plot!(positive_freqs, positive_magnitude, xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum",title="Fourier Transform", yscale=:log10)
    #plot!(freqs, abs.(fft_result), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum",color=:red,xlim=(0,0.00045))
end
plot!(legend=:false)
vline!([0.000085],color=:black,legend=:false)

plot();
for i in 1:18
    fft_result = fft(Float64.(OG[:,i]))
    magnitude = abs.(fft_result) ./ 216
    freqs = fftfreq(length(time),0.00083333)
    positive_indices = findall(>=(0), freqs)
    positive_freqs = freqs[positive_indices]
    positive_magnitude = magnitude[positive_indices]
    plot!(positive_freqs, positive_magnitude, xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum", yscale=:log10)
    #plot!(freqs, abs.(fft_result), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum",color=:red,xlim=(0,0.00045))
end
plot!(legend=:false)
vline!([0.000085],color=:black,legend=:false)

per=[];
freqsp=[];
#PSD
plot();
for i in 1:18
    p = periodogram(Float64.(detrended_datam[:,i]), fs=sr)
    freqsp = p.freq
    psd = p.power
    push!(freqsp_d,freqsp)
    push!(per,psd)
    positive_indices = findall(>(0), freqsp)
    positive_freqsp = freqsp[positive_indices]
    positive_psd = psd[positive_indices]
    plot!(positive_freqsp,positive_psd, xlabel="Frequency (Hz)", ylabel="Power Spectral Density",
    title="Power Spectral Density", legend=false, yscale=:log10)
end
plot!(legend=:false)

plot();
max_lag = 215 
lags = 0:215 # Maximum number of lags to compute
xl = collect(0:max_lag)
for i in 1:18
    autocorrelation = autocor(Float64.(detrended_datam[:,i]), xl)
    plot!(xl, autocorrelation)
end
plot!(legend=:false,xlabel="Lag", ylabel="Autocorrelation", title="Autocorrelation (Detrended + Mean)",ylims=(-0.3,1))


plot();
max_lag = 215 
lags = 0:215 # Maximum number of lags to compute
xl = collect(0:max_lag)
for i in 1:18
    autocorrelation = autocor(Float64.(OG[:,i]), xl)
    plot!(xl, autocorrelation)
end
plot!(legend=:false,xlabel="Lag", ylabel="Autocorrelation", title="Autocorrelation Raw Data")



fft_result = fft(OG1)
freqs = fftfreq(length(time), 0.00083333)
plot(freqs[1:div(end,2)], abs.(fft_result[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum")
plot();
for i in 1:18
    fft_result = fft(Float64.(OG[:,i]))
    freqs = fftfreq(length(time), 0.00083333)
    plot!(freqs[1:div(end,2)], abs.(fft_result[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum",color=:red)
end
plot!(legend=:false, xlimit=(-0.00001,0.0001))
for i in 1:18
    fft_result = fft(Float64.(Cor[:,i]))
    freqs = fftfreq(length(time), 0.00083333)
    plot!(freqs[1:div(end,2)], abs.(fft_result[1:div(end,2)]), xlabel="Frequency (Hz)", ylabel="Magnitude", label="Frequency Spectrum", color=:black)
end
plot!(legend=:false, xlimit=(-0.00001,0.0001))
vline!([0.000005],color=:black,legend=:false)


lp = Lowpass(0.000085, fs=0.00083333)
dmeth = Butterworth(2)
OG_means = [];
filt_OG_d3 = [];
filt_OG_d85 = [];
fluct_OG_d = [];
for i in 1:18
    k1 = mean(Float64.(detrended_datam[:,i]))
    push!(OG_means,k1)
    filt_OG = filtfilt(digitalfilter(lp,dmeth), Float64.(detrended_datam[:,i]))
    push!(filt_OG_d85,filt_OG)
    fluct_OG = Float64.(detrended_datam[:,i]) - filt_OG
    push!(fluct_OG_d,fluct_OG)
end
filt_OG_traces85 = Float64.(DataFrame([Symbol("Col$i") => vec for (i, vec) in enumerate(filt_OG_d85)]))


plot(time,fluct_OG_d,legend=:false)
plot();
for i in 1:Cols
    plot!(time,filt_OG_traces3[:,i]/1000000,color=:black,legend=:false,title=:"Low Pass Filtered Traces\nand Linear Detrended + Mean Traces", xlabel="Time (Hrs)", ylabel = "# GFP Molecules (10^6)")
    plot!(time,detrended_datam[:,i]/1000000,legend=:false)
end
plot!()

plot!(time,filt_OG_traces3[:,1]/1000000,color=:black,linestyle=:dash,legend=:false,title=:"Low Pass filtered traces\nand linear detrended + Mean traces", xlabel="Time (Hrs)", ylabel = "# GFP Molecules (10^6)")
plot!(time,filt_OG_traces85[:,1]/1000000,color=:black,legend=:false,title=:"Low Pass filtered traces\nand linear detrended + Mean traces", xlabel="Time (Hrs)", ylabel = "# GFP Molecules (10^6)")
plot!(time,detrended_datam[:,1]/1000000,color=:red,legend=:false)
plot();
#Get Experimental time trace data High Expression Example
BFCexample_GFP = XLSX.readxlsx("C:/Users/Strey Lab/Documents/GitHub/Analysis_of_SGC_Simulations/StreyCats/Broken Feedback Circuit/Broken_Feedback_Circuit_High_Expression_Example.xlsx");
GFPexample = BFCexample_GFP["Sheet1"];
GFPexample = GFPexample[:];
example_GFP = GFPexample[2:end,2:end]/1000000

plot!(time, example_GFP)
plot!(time, OG[:,1]/1000000)
plot();

fft_resultB = fft(Float64.(detrended_datam[:,1]))
freqsB = fftfreq(length(time), 0.00083333)
magnitudeB = abs.(fft_resultB) ./ 216
fft_resultE = fft(example_GFP*1000000)
freqsE = fftfreq(length(time), 0.00083333)
magnitudeE = abs.(fft_resultE) ./ 216
positive_indicesB = findall(>=(0), freqsB)
positive_freqsB = freqsB[positive_indicesB]
positive_fftB = magnitudeB[positive_indicesB]
positive_fftBN =positive_fftB/positive_fftB[1]
positive_indicesE = findall(>=(0), freqsE)
positive_freqsE = freqsE[positive_indicesE]
positive_fftE = magnitudeE[positive_indicesE]
positive_fftEN = positive_fftE/positive_fftE[1]

plot!(positive_freqsB, positive_fftBN, xlabel="Frequency (Hz)", ylabel="Magnitude", label="BKG Frequency Spectrum",color=:red, yscale=:log10)
plot!(positive_freqsE, positive_fftEN, xlabel="Frequency (Hz)", ylabel="Magnitude", label="Cell Frequency Spectrum",color=:blue, yscale=:log10 )
vline!([0.00005],color=:black,legend=:false)
plot();

dif_fft = magnitudeE - magnitudeB
positive_indicesD = findall(>(0), freqsB)
positive_freqsD = freqsB[positive_indicesD]
positive_fftD = dif_fft[positive_indicesD]
plot(positive_freqsD, positive_fftD, xlabel="Frequency (Hz)", ylabel="Magnitude", label="BKG Frequency Spectrum",color=:black)

plot();
dif_t = real(ifft(dif_fft))
plot!(time, example_GFP*1000000)
plot!(time, dif_t)


#Fourier Transform of mean of experimental data and mean of background data
BFC_GFP = XLSX.readxlsx("C:/Users/Strey Lab/Documents/GitHub/Analysis_of_SGC_Simulations/StreyCats/Broken Feedback Circuit/Broken Feedback Circuit Time Traces 5_21_24 GFP.xlsx");
GFPex = BFC_GFP["Sheet1"];
GFPex = GFPex[:];
ex_GFP = GFPex[2:end,2:end];
ex_GFP_m = mean(ex_GFP, dims=2)
plot(time, ex_GFP_m)
plot();
for i in 1:965
    plot!(time,ex_GFP[:,i]/1000000)
end
plot!(time, ex_GFP_m/1000000, color=:black, legend=:false,xlabel="Time (hrs)",ylabel="#GFP (x10^6)")

BKG_OG_det_m = mean(detrended_datam, dims=2)
filt_BKG_OG_det_m = filtfilt(digitalfilter(lp,dmeth), Float64.(BKG_OG_det_m))
plot!(time, BKG_OG_det_m/1000000,xlabel="Time (hrs)",ylabel="#GFP (x10^6)", color=:red, label="Mean of Background Traces")
plot!(time, filt_BKG_OG_det_m/1000000,xlabel="Time (hrs)",ylabel="#GFP (x10^6)", color=:Black, label=:"LP: 8.5x10^-5 Hz")
plot();
for i in 1:18
    plot!(time, detrended_datam[:,i]/1000000)
end
plot!(time, BKG_OG_det_m/1000000, color=:black, legend=:false,xlabel="Time (hrs)",ylabel="#GFP (x10^6)")
plot();
fft_resultB = fft(Float64.(BKG_OG_det_m))
freqsB = fftfreq(length(time), 0.00083333)
magnitudeB = abs.(fft_resultB) ./ 216
fft_resultE = fft(ex_GFP_m)
freqsE = fftfreq(length(time), 0.00083333)
magnitudeE = abs.(fft_resultE) ./ 216
positive_indicesB = findall(>=(0), freqsB)
positive_freqsB = freqsB[positive_indicesB]
positive_fftB = magnitudeB[positive_indicesB]
positive_fftBN =positive_fftB/positive_fftB[1]
positive_indicesE = findall(>=(0), freqsE)
positive_freqsE = freqsE[positive_indicesE]
positive_fftE = magnitudeE[positive_indicesE]
positive_fftEN = positive_fftE/positive_fftE[1]

plot!(positive_freqsB, positive_fftBN, xlabel="Frequency (Hz)", ylabel="Magnitude", label="Mean BKG Frequency Spectrum",color=:red, yscale=:log10)
plot!(positive_freqsE, positive_fftEN, xlabel="Frequency (Hz)", ylabel="Magnitude", label="Mean Cell Frequency Spectrum",color=:blue, yscale=:log10 )
vline!([0.000085],color=:black, label="LP: 8.5x10^-5 Hz")
vline!([0.0000333], color=:black, linestyle=:dash, label="LP: 3.3x10^-5 Hz")
plot();