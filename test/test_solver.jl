using Test
using EtFsp
using LinearAlgebra
using Sundials: CVODE_BDF

## Bursting gene model 
𝕊 = [[-1, 1, 0] [1, -1, 0] [0, 0, 1] [0, 0, -1]]
x₀ = [1, 0, 0]
k₀₁ = 0.05
k₁₀ = 0.1
λ = 5.0
γ = 1.0

propensities = [
    propensity((G₀::Int, G₁::Int, RNA::Int) -> k₀₁ * G₀)
    propensity_timevarying(t -> max(0.0, 1.0 - sin(π * t / 2)), (G₀::Int, G₁::Int, RNA::Int) -> k₁₀ * G₁)
    propensity((G₀::Int, G₁::Int, RNA::Int) -> λ * G₁)
    propensity((G₀::Int, G₁::Int, RNA::Int) -> γ * RNA)
]
bursting_model = CmeModel(𝕊, propensities)

propensities_parmaterized = [
    propensity((G₀::Int, G₁::Int, RNA::Int, k₀₁, k₁₀, λ, γ) -> k₀₁ * G₀)
    propensity_timevarying((t, k₀₁, k₁₀, λ, γ) -> max(0.0, 1.0 - sin(π * t / 2)), (G₀::Int, G₁::Int, RNA::Int,  k₀₁, k₁₀, λ, γ) -> k₁₀ * G₁)
    propensity((G₀::Int, G₁::Int, RNA::Int, k₀₁, k₁₀, λ, γ) -> λ * G₁)
    propensity((G₀::Int, G₁::Int, RNA::Int, k₀₁, k₁₀, λ, γ) -> γ * RNA)
]
bursting_model_parameterized = CmeModel(𝕊, propensities_parmaterized)

# Fixed FSP solver 
# Solve to get dense outputs
tspan = (0.0, 120.0)
𝔛 = SparseStateSpace(𝕊, x₀)
expand!(𝔛, 20)
p0 = SparseMultIdxVector(𝔛, [[1, 0, 0] => 1.0])
fspmethod = FixedSparseFsp(CVODE_BDF(linear_solver = :GMRES))
solutions = solve(bursting_model, p0, tspan, fspmethod, odertol = 1.0e-4, odeatol = 1.0e-14);
@test prod(
    [(sum(p) + sum(sinks) ≈ 1.0) for (p, sinks) in zip(solutions.p, solutions.sinks)
])

# Solve to get outputs at specific times
tspan = (0.0, 120.0)
toutputs = 0.0:20.0:120.0
𝔛 = SparseStateSpace(𝕊, x₀)
expand!(𝔛, 20)
p0 = SparseMultIdxVector(𝔛, [[1, 0, 0] => 1.0])
fspmethod = FixedSparseFsp(CVODE_BDF(linear_solver = :GMRES))
solutions = solve(bursting_model, p0, tspan, fspmethod, odertol = 1.0e-4, odeatol = 1.0e-14, saveat = toutputs);
@test prod(
    [(sum(p) + sum(sinks) ≈ 1.0) for (p, sinks) in zip(solutions.p, solutions.sinks)
])
@test length(solutions) == length(toutputs)

# Consistency between using parameter-free representation and parametric representation
tspan = (0.0, 120.0)
toutputs = 0.0:20.0:120.0
p0 = SparseMultIdxVector(SparseStateSpace(𝕊, x₀), [x₀ => 1.0])
fspmethod = AdaptiveSparseFsp(
    ode_method = CVODE_BDF(linear_solver = :GMRES),
    space_adapter = SelectiveRStepAdapter(10, 10)
)

fspsolutions1 = solve(bursting_model, p0, tspan, fspmethod, odertol = 1.0e-4, odeatol = 1.0e-14, saveat = toutputs)
@test prod(
    [(sum(p) + sum(sinks) ≈ 1.0) for (p, sinks) in zip(fspsolutions1.p, fspsolutions1.sinks)
])

fspsolutions2 = solve(bursting_model_parameterized, p0, tspan, fspmethod, [k₀₁, k₁₀, λ, γ], odertol = 1.0e-4, odeatol = 1.0e-14, saveat = toutputs)
@test prod(
    [(sum(p) + sum(sinks) ≈ 1.0) for (p, sinks) in zip(fspsolutions2.p, fspsolutions2.sinks)
])

@test prod(
    [norm(p1.values - p2.values, 1) ≤ 1.0E-14 for (p1, p2) in zip(fspsolutions1.p, fspsolutions2.p)]
)