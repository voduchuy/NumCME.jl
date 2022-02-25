using Test 
using EtFsp 
using Sundials: CVODE_BDF 

## Bursting gene model 
𝕊 = [[-1, 1, 0] [1, -1, 0] [0, 0, 1] [0, 0, -1]]
x₀ = [1, 0, 0]
k₀₁ = 0.05
k₁₀ = 0.1
λ = 5.0
γ = 1.0

# TimeSeparablePropensity functions for the time-varying case 
propensities_tv = [
    TimeSeparablePropensity(nothing, (G₀::Int, G₁::Int, RNA::Int) -> k₀₁ * G₀)
    TimeSeparablePropensity(t -> max(0.0, 1.0 - sin(π * t / 2)), (G₀::Int, G₁::Int, RNA::Int) -> k₁₀ * G₁)
    TimeSeparablePropensity(nothing, (G₀::Int, G₁::Int, RNA::Int) -> λ * G₁)
    TimeSeparablePropensity(nothing, (G₀::Int, G₁::Int, RNA::Int) -> γ * RNA)
]
bursting_model_tvarying = CmeModel(𝕊, propensities_tv)

# Fixed FSP solver 
# Solve to get dense outputs
tspan = (0.0, 120.0)
𝔛 = SparseStateSpace(𝕊, x₀)
expand!(𝔛, 20)
p0 = FspSparseVector{Float64}(𝔛, [[1,0,0]=>1.0])
fspmethod = FixedSparseFsp(CVODE_BDF(linear_solver=:GMRES))
solutions = solve(bursting_model_tvarying, p0, tspan, fspmethod, odertol=1.0e-4, odeatol=1.0e-14);
@test prod([(sum(sol.p) + sum(sol.sinks) ≈ 1.0) for sol in solutions])

# Solve to get outputs at specific times
tspan = (0.0, 120.0)
toutputs = 0.0:20.0:120.0
𝔛 = SparseStateSpace(𝕊, x₀)
expand!(𝔛, 20)
p0 = FspSparseVector{Float64}(𝔛, [[1,0,0]=>1.0])
fspmethod = FixedSparseFsp(CVODE_BDF(linear_solver=:GMRES))
solutions = solve(bursting_model_tvarying, p0, tspan, fspmethod, odertol=1.0e-4, odeatol=1.0e-14, saveat=toutputs);
@test prod([(sum(sol.p) + sum(sol.sinks) ≈ 1.0) for sol in solutions])
@test length(solutions) == length(toutputs)

