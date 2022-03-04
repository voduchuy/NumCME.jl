using EtFsp
using BenchmarkTools
import DifferentialEquations as DE
using Sundials: CVODE_BDF

##  Toggle-switch model
𝕊 = [[1, 0] [-1, 0] [0, 1] [0, -1]]
α₁ = propensity((S₁, S₂, b₁, b₂, k₁, k₂, a₂₁, a₁₂, n₂₁, n₁₂, γ₁, γ₂, UV, Δtᵤᵥ) -> b₁ + k₁/(1.0 + a₂₁*S₂^n₂₁))
α₂ = propensity((S₁, S₂, b₁, b₂, k₁, k₂, a₂₁, a₁₂, n₂₁, n₁₂, γ₁, γ₂, UV, Δtᵤᵥ) -> γ₁ * S₁)
α₃ = propensity((S₁, S₂, b₁, b₂, k₁, k₂, a₂₁, a₁₂, n₂₁, n₁₂, γ₁, γ₂, UV, Δtᵤᵥ) -> b₂ + k₂/(1.0 + a₁₂*S₁^n₁₂))

# Mathematically equivalent definitions of the fourth propensity function for S2 degradation, but computationally they are different: α₄ requires the time and state variables to be evaluated jointly whereas β₄ is factored into a time-only and a state-only functions
α₄ = propensity_timevarying((t, S₁, S₂, b₁, b₂, k₁, k₂, a₂₁, a₁₂, n₂₁, n₁₂, γ₁, γ₂, UV, Δtᵤᵥ) -> (γ₂ + (t ≤ Δtᵤᵥ)*0.002*UV^2/(1260+UV^3))*S₂)
β₄ = propensity_timevarying((t, b₁, b₂, k₁, k₂, a₂₁, a₁₂, n₂₁, n₁₂, γ₁, γ₂, UV, Δtᵤᵥ) -> γ₂ + (t ≤ Δtᵤᵥ)*0.002*UV^2/(1260+UV^3), (S₁, S₂, b₁, b₂, k₁, k₂, a₂₁, a₁₂, n₂₁, n₁₂, γ₁, γ₂, UV, Δtᵤᵥ) -> S₂)

propensities_joint = [α₁,α₂,α₃,α₄]
propensities_separable = [α₁,α₂,α₃,β₄]
model_joint = CmeModel(𝕊, propensities_joint)
model_separable = CmeModel(𝕊, propensities_separable)

x₀ = [0, 0]

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

𝔛₀ = SparseStateSpace(𝕊, x₀)
p0 = SparseMultIdxVector(𝔛₀, [x₀=>1.0])

tspan = (0.0, 8.0*3600)
saveat = 0.0:60.0:8*3600.0

fixedrstepfsp = AdaptiveSparseFsp(
    ode_method = CVODE_BDF(linear_solver=:GMRES),
    space_adapter = RStepAdapter(20, 5, true)
)

adaptiverstepfsp = AdaptiveSparseFsp(
    ode_method = CVODE_BDF(linear_solver=:GMRES),
    space_adapter = SelectiveRStepAdapter(20, 5, true)
)

@btime fspsol1 = solve(model_separable, p0, tspan, fixedrstepfsp, θ, saveat=saveat, odertol=1.0E-4, odeatol=1.0E-14);
@btime fspsol2 = solve(model_separable, p0, tspan, adaptiverstepfsp, θ, saveat=saveat, odertol=1.0E-4, odeatol=1.0E-14);
@btime fspsol3 = solve(model_joint, p0, tspan, fixedrstepfsp, θ, saveat=saveat, odertol=1.0E-4, odeatol=1.0E-14);
@btime fspsol4 = solve(model_joint, p0, tspan, adaptiverstepfsp, θ, saveat=saveat, odertol=1.0E-4, odeatol=1.0E-14);

