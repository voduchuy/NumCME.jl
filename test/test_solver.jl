using Test
using Julifsp
using LinearAlgebra
using Sundials: CVODE_BDF

## Bursting gene model 
𝕊 = [[-1, 1, 0] [1, -1, 0] [0, 0, 1] [0, 0, -1]]
x₀ = [1, 0, 0]
k₀₁ = 0.05
k₁₀ = 0.1
λ = 5.0
γ = 1.0

θ = [k₀₁, k₁₀, λ, γ]

# Propensity formulations that have no dependencies on parameters
a1 = propensity() do x, p
    k₀₁ * x[1]
end
a2 = propensity_timevarying((t, p) -> max(0.0, 1.0 - sin(π * t / 2))) do x, p
    k₁₀ * x[2]
end
a2j = propensity_timevarying() do t, x, p
    max(0.0, 1.0 - sin(π * t / 2)) * k₁₀ * x[2]
end
a3 = propensity() do x, p
    λ * x[2]
end
a4 = propensity() do x, p
    γ * x[3]
end

# Propensity formulations that have dependencies on parameters
a1_p = propensity() do x, p
    p[1] * x[1]
end
a2_p = propensity_timevarying((t, p) -> max(0.0, 1.0 - sin(π * t / 2))) do x, p
    p[2] * x[2]
end
a3_p = propensity() do x, p
    p[3] * x[2]
end
a4_p = propensity() do x, p
    p[4] * x[3]
end


bursting_model = CmeModel(𝕊, [a1, a2, a3, a4], [])
bursting_model_parameterized = CmeModel(𝕊, [a1_p, a2_p, a3_p, a4_p], θ)

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
    space_adapter = SelectiveRStepAdapter(10, 10, true)
)

fspsolutions1 = solve(bursting_model, p0, tspan, fspmethod, odertol = 1.0e-4, odeatol = 1.0e-14, saveat = toutputs)
@test prod(
    [(sum(p) + sum(sinks) ≈ 1.0) for (p, sinks) in zip(fspsolutions1.p, fspsolutions1.sinks)
])

fspsolutions2 = solve(bursting_model_parameterized, p0, tspan, fspmethod, odertol = 1.0e-4, odeatol = 1.0e-14, saveat = toutputs)
@test prod(
    [(sum(p) + sum(sinks) ≈ 1.0) for (p, sinks) in zip(fspsolutions2.p, fspsolutions2.sinks)
])

@test prod(
    [norm(p1.values - p2.values, 1) ≤ 1.0E-14 for (p1, p2) in zip(fspsolutions1.p, fspsolutions2.p)]
)