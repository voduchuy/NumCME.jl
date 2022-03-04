using Julifsp
using Test
using LinearAlgebra: norm 


## Bursting gene model 
𝕊 = [[-1, 1, 0] [1, -1, 0] [0, 0, 1] [0, 0, -1]]
x₀ = [1, 0, 0]
k₀₁ = 0.05
k₁₀ = 0.1
λ = 5.0
γ = 1.0
# TimeSeparablePropensity functions for the time-invariant case
propensities_tind = [
    propensity((G₀::Int, G₁::Int, RNA::Int) -> k₀₁ * G₀)
    propensity((G₀::Int, G₁::Int, RNA::Int) -> k₁₀ * G₁)
    propensity((G₀::Int, G₁::Int, RNA::Int) -> λ * G₁)
    propensity((G₀::Int, G₁::Int, RNA::Int) -> γ * RNA)
]

# TimeSeparablePropensity functions for the time-varying case 
propensities_tv = [
    propensity((G₀::Int, G₁::Int, RNA::Int) -> k₀₁ * G₀)
    propensity_timevarying(t -> max(0.0, 1.0 - sin(π * t / 2)), (G₀::Int, G₁::Int, RNA::Int) -> k₁₀ * G₁)
    propensity((G₀::Int, G₁::Int, RNA::Int) -> λ * G₁)
    propensity((G₀::Int, G₁::Int, RNA::Int) -> γ * RNA)
]

bursting_model_tinvar = CmeModel(𝕊, propensities_tind)
bursting_model_tvarying = CmeModel(𝕊, propensities_tv)

𝔛 = SparseStateSpace(𝕊, x₀)
expand!(𝔛, 2)
𝐀 = FspMatrixSparse(𝔛, propensities_tind)
@test size(𝐀, 1) == get_state_count(𝔛) + get_sink_count(𝔛)
@test size(𝐀, 2) == get_state_count(𝔛) + get_sink_count(𝔛)
𝐯 = ones(Float64, size(𝐀, 1))
𝐰 = 𝐀(1.0) * 𝐯
@test sum(𝐰) ≈ 0.0 atol = 1.0e-14
𝐰 = 𝐀 * 𝐯
@test sum(𝐰) ≈ 0.0 atol = 1.0e-14

# Test mat-vec for time-varying matrix
𝔛 = SparseStateSpace(𝕊, x₀)
expand!(𝔛, 2)
A1 = FspMatrixSparse(𝔛, propensities_tv)
@test size(A1, 1) == get_state_count(𝔛) + get_sink_count(𝔛)
@test size(A1, 2) == get_state_count(𝔛) + get_sink_count(𝔛)
𝐯 = ones(Float64, size(A1, 1))
w1 = A1(1.0) * 𝐯
@test sum(w1) ≈ 0.0 atol = 1.0e-14
w1 = A1 * 𝐯
@test sum(w1) ≈ 0.0 atol = 1.0e-14

propensities_tv2 = [
    propensity((G₀::Int, G₁::Int, RNA::Int) -> k₀₁ * G₀)
    propensity_timevarying((t, G₀::Int, G₁::Int, RNA::Int) -> max(0.0, 1.0 - sin(π * t / 2))* k₁₀ * G₁)
    propensity((G₀::Int, G₁::Int, RNA::Int) -> λ * G₁)
    propensity((G₀::Int, G₁::Int, RNA::Int) -> γ * RNA)
]
A2 = FspMatrixSparse(𝔛, propensities_tv2)
w2 = A2(1.0) * 𝐯
@test sum(𝐰) ≈ 0.0 atol = 1.0e-14
w2 = A2 * 𝐯
@test sum(w2) ≈ 0.0 atol = 1.0e-14

@test norm(w1 -w2) ≈ 0
# ## Test the integration of the FSP system 
# using DifferentialEquations: ODEProblem
# import Sundials 
# using SparseArrays
# 𝔛 = SparseStateSpace(𝕊, x₀)
# expand!(𝔛, 20)
# 𝐀 = FspMatrixSparse(𝔛, propensities_tv)
# 𝐩₀ = [1.0;zeros(Float64, size(𝐀, 1) - 1)]
# tspan = (0.0, 120.0)
# function fsprhs!(du, u, θ, t)
#     matvec!(du, t, 𝐀, u)    
#     nothing 
# end
# fspprob = ODEProblem(fsprhs!, 𝐩₀, tspan)
# @time sol = Sundials.solve(fspprob, Sundials.CVODE_BDF(linear_solver=:GMRES), atol=1.0e-14, rtol=1.0e-4);
# @test prod([sum(p) ≈ 1.0  for p in sol.u])


