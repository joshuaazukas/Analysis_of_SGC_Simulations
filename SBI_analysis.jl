using JSON3
using Plots
json_string = read("real_results_bfc_dM.json", String)

hello_world = JSON3.read(json_string)

m = hello_world["means"]
s = hello_world["standard_deviations"]
mm = reduce(hcat,m)
ss = reduce(hcat,s)

histogram(mm[1,:],label = "kon",bins=50)
savefig("kon.png")
histogram(mm[2,:],label = "koff",bins=50)
savefig("koff.png")
histogram(mm[3,:],label = "kTr",bins=50)
savefig("kTr.png")
histogram(mm[4,:],label = "KTl",bins=50)
savefig("kTl.png")
histogram(mm[5,:],label = "dG",bins=50)
savefig("dG.png")
histogram(mm[6,:],label = "dM",bins=50)
savefig("dM.png")
scatter(mm[3,:],mm[4,:],xlim=(0,900), label="kTr vs kTl")
savefig("kTrKTl.png")

scatter(mm[1,:],mm[2,:], label="kon vs koff")
savefig("konKoff.png")

scatter(mm[3,:],mm[4,:], label="kTr vs kTl")
savefig("kTrKTl.png")

scatter(mm[5,:],mm[6,:], label="dG vs dM")
savefig("dGdM.png")