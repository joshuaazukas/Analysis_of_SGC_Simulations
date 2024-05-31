using Catalyst, DifferentialEquations, Plots, DataFrames, CSV, Interpolations
plt1 = []
plt2 = []
plt3 = []
plt4 = []
plt5 = []
plt6 = []
plt7 = []
t_inter = 0:0.1:300
for i in 1:4 # number of time the loop will run
    rn = @reaction_network begin
        (kOn,kOff), A + DNA <--> A_DNA  # Activator binding DNA
        (kOnt,kOfft), A_DNA + GFP <--> A_DNA_T # TetR/GFP Binds active promoter, prevent transcription
        (kOnt,kOfft), DNA + GFP <--> DNA_T  # TetR/GFP Binds inactive promoter, prevent transcription
        (kOn,kOff), A + DNA_T <--> A_DNA_T # Activator binds TetR bound promoter, prevent transcription
        k, A_DNA --> A_DNA + RNA # Transcription of DNA with Activator Bound
        kT, RNA --> RNA + GFP # Translation of RNA to produce TetR and GFP
        deg_R, RNA --> 0 # Degradation of RNA
        deg_G, GFP --> 0 # Degradation of protein
    end kOn kOff kOnt kOfft k kT deg_R deg_G

    #:kOnt => 0.00000001*60*60, :kOfft => 0.0002*60*60,

    tspan = (0.0, 300)
    u0 = [:A => 1, :DNA => 1, :A_DNA => 0, :DNA_T => 0, :A_DNA_T => 0, :RNA => 0, :GFP => 0]
    p = [:kOn => 0.0042*60*60, :kOff => 0.0038*60*60, :kOnt => 0, :kOfft => 0, :k => 250, :kT => 10.0, :deg_R => 0.8, :deg_G => 0.3]

    dprob = DiscreteProblem(rn, u0, tspan, p)
    jprob = JumpProblem(rn, dprob, Direct())
    @time sol = solve(jprob, SSAStepper())

    Plots.plot(sol)
    Plots.plot(sol,idxs=1, xlimit=(0.0,1.0))
    Plots.plot(sol,idxs=2, xlimit=(0.0,1.0))
    Plots.plot(sol,idxs=3, xlimit=(0.0,1.0))
    Plots.plot(sol,idxs=4)
    Plots.plot(sol,idxs=5)

    a = (sol.t, sol[1,:])
    append!(plt1,a)

    dna = (sol.t, sol[2,:])
    push!(plt2,dna)

    a_dna = (sol.t, sol[3,:])
    push!(plt3,a_dna)

    a_dna_t =  (sol.t, sol[5,:])
    push!(plt5,a_dna_t)

    dna_t = (sol.t, sol[6,:])
    push!(plt4,dna_t)

    rna_int = linear_interpolation(sol.t, sol[7,:])
    rna = rna_int.(t_inter)
    push!(plt6,rna)
    
    gfp_int = linear_interpolation(sol.t, sol[4,:])
    gfp = gfp_int.(t_inter)
    push!(plt7,gfp)
end

Plots.plot(plt7)