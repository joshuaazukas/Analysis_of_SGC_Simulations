using DifferentialEquations
using ModelingToolkit
using Plots

@variables t

D = Differential(t)

@parameters begin
    k_OU
    b_OU
    kT
    deg_R
    deg_G
    μ
end

@variables begin
    k_RNA
    rna
    gfp
end

@brownian w

# For the OU process <xx(t)> = b_OU^2/(2*k_OU)* exp(-k_OU*t)
eqs =  [
    D(k_RNA) ~ -k_OU * (k_RNA - μ)  + b_OU * w
    D(rna) ~ k_RNA - deg_R * rna
    D(gfp) ~ rna - deg_G * gfp ]

# to calculate the parameters we need to know the mean and variance of the OU process
# resulting from master equation describing the activator complex binding
# λ 0->1 ω 1->0 transition rates 
# μ = k * λ/(λ + ω)
# var = k^2 ω*λ/(λ + ω)^2
# b_OU = (λ + ω)

p = []
using DifferentialEquations: solve
@mtkbuild rnasde = RNASDE()

# Create the system
sde_system = SDEProblem(eqs, [k_RNA => 0.0, rna => 0.0, gfp => 0.0], t, [k_OU => 1.0, b_OU => 0.1, deg_R => 0.1, deg_G => 0.1, μ => 5.0], noise=[w])

# Solve the system
sol = solve(sde_system, EM())

# Plot the solution
plot(sol, vars=[k_RNA, rna, gfp], xlabel="Time", ylabel="Concentration", legend=:topright)