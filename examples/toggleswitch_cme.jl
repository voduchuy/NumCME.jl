using NumCME
using BenchmarkTools
import DifferentialEquations as DE
using Sundials: CVODE_BDF
using LinearAlgebra: BLAS

##  Toggle-switch model
𝕊 = [[1, 0] [-1, 0] [0, 1] [0, -1]]

α₁ = propensity() do x, p
    S₁, S₂ = x
    b₁, b₂, k₁, k₂, a₂₁, a₁₂, n₂₁, n₁₂, γ₁, γ₂, UV, Δtᵤᵥ = p
    b₁ + k₁ / (1.0 + a₂₁ * S₂^n₂₁)
end
α₂ = propensity() do x, p
    S₁, S₂ = x
    b₁, b₂, k₁, k₂, a₂₁, a₁₂, n₂₁, n₁₂, γ₁, γ₂, UV, Δtᵤᵥ = p
    γ₁ * S₁
end
α₃ = propensity() do x, p
    S₁, S₂ = x
    b₁, b₂, k₁, k₂, a₂₁, a₁₂, n₂₁, n₁₂, γ₁, γ₂, UV, Δtᵤᵥ = p
    b₂ + k₂ / (1.0 + a₁₂ * S₁^n₁₂)
end
α₄ = propensity() do t, x, p
    S₁, S₂ = x
    b₁, b₂, k₁, k₂, a₂₁, a₁₂, n₂₁, n₁₂, γ₁, γ₂, UV, Δtᵤᵥ = p
    (γ₂ + (t ≤ Δtᵤᵥ) * 0.002 * UV^2 / (1260 + UV^3)) * S₂
end

# This propensity formulation is mathematically equivalent to α₄ but leads to more computationally efficient CME solves because β₄ is factored into a time-only and a state-only functions
function degradation_rate(t, p)
    b₁, b₂, k₁, k₂, a₂₁, a₁₂, n₂₁, n₁₂, γ₁, γ₂, UV, Δtᵤᵥ = p
    γ₂ + (t ≤ Δtᵤᵥ) * 0.002 * UV^2 / (1260 + UV^3)
end
β₄ = propensity(degradation_rate) do x, p
    x[2]
end

b₁ = 2.2E-3
b₂ = 6.8E-5
k₁ = 1.7E-2
k₂ = 1.6E-2
a₂₁ = 2.6E-3
a₁₂ = 6.1E-3
n₂₁ = 3
n₁₂ = 2.1
γ₁ = 3.8E-4
γ₂ = 3.8E-4
UV = 10.0
Δtᵤᵥ = 3600

θ = [b₁, b₂, k₁, k₂, a₂₁, a₁₂, n₂₁, n₁₂, γ₁, γ₂, UV, Δtᵤᵥ]

propensities_joint = [α₁, α₂, α₃, α₄]
propensities_separable = [α₁, α₂, α₃, β₄]
model_joint = CmeModel(𝕊, propensities_joint, θ)
model_separable = CmeModel(𝕊, propensities_separable, θ)

x₀ = [0, 0]
𝔛₀ = StateSpaceSparse(𝕊, x₀)
p0 = FspVectorSparse(𝔛₀, [x₀ => 1.0])

tspan = (0.0, 8.0 * 3600)
saveat = 0.0:60.0:8*3600.0

fixedrstepfsp = AdaptiveFspSparse(
    ode_method = CVODE_BDF(linear_solver = :GMRES),
    space_adapter = RStepAdapter(20, 5, true)
)

adaptiverstepfsp = AdaptiveFspSparse(
    ode_method = CVODE_BDF(linear_solver = :GMRES),
    space_adapter = SelectiveRStepAdapter(20, 5, true)
)

for threadcount in [1, 4, 8]
    BLAS.set_num_threads(threadcount)
    println("Solving with full R-step expansion and separable propensity format")
    @btime fspsol1 = solve(model_separable, p0, tspan, fixedrstepfsp, saveat = saveat, odertol = 1.0E-4, odeatol = 1.0E-14)
    println("Solving with selective R-step expansion and separable propensity format")
    @btime fspsol2 = solve(model_separable, p0, tspan, adaptiverstepfsp, saveat = saveat, odertol = 1.0E-4, odeatol = 1.0E-14)
    println("Solving with full R-step expansion and non-separable propensity format")
    @btime fspsol3 = solve(model_joint, p0, tspan, fixedrstepfsp, saveat = saveat, odertol = 1.0E-4, odeatol = 1.0E-14)
    println("Solving with selective R-step expansion and non-separable propensity format")
    @btime fspsol4 = solve(model_joint, p0, tspan, adaptiverstepfsp, saveat = saveat, odertol = 1.0E-4, odeatol = 1.0E-14)
end



