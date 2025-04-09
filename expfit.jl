using LsqFit, Plots, XLSX

data_reg = XLSX.readxlsx("LightInd_TL2Reg_Time_Traces.xlsx")
sh1 = data_reg["Sheet1"]

traces = Float64.(sh1["A2:ALB217"])

t_length = size(traces)[1]
cell_length = size(traces)[2]

t = collect(1:t_length)/3.0

plot(t,traces[:,:], legend=nothing, xlabel="hours")

expfit(x,p) = p[1] .* exp.(p[2] .* x)

res = curve_fit(expfit,t,traces[:,1],[1e6,1/60])
res.param

param_list = []
for i in 1:cell_length
    res = curve_fit(expfit,t,traces[:,i],[1e6,1/60])
    push!(param_list,res.param)
end

params = reduce(hcat,param_list)'

histogram(params[:,2])
savefig("expdist.png")