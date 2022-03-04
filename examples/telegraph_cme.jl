using EtFsp
import DifferentialEquations as DE
using Sundials: CVODE_BDF
using StaticArrays

## Bursting gene model 
𝕊 = [[-1, 1, 0] [1, -1, 0] [0, 0, 1] [0, 0, -1]]
x₀ = [1, 0, 0]
k₀₁ = 0.05
k₁₀ = 0.1
λ = 5.0
γ = 0.5

propensities_tv = [
    propensity((G₀, G₁, RNA, k₀₁, k₁₀, λ, γ) -> k₀₁ * G₀)
    propensity((G₀, G₁, RNA, k₀₁, k₁₀, λ, γ) -> k₁₀ * G₁)
    propensity((G₀, G₁, RNA, k₀₁, k₁₀, λ, γ) -> λ * G₁)
    propensity((G₀, G₁, RNA, k₀₁, k₁₀, λ, γ) -> γ * RNA)
]

θ = [k₀₁, k₁₀, λ, γ]
model = CmeModel(𝕊, propensities_tv)
𝔛₀ = SparseStateSpace(model.stoich_matrix, x₀)
expand!(𝔛₀, 10)
p0 = SparseMultIdxVector(𝔛₀, [x₀=>1.0])
tspan = (0.0, 300.0)

fixedrstepfsp = AdaptiveSparseFsp(
    ode_method = CVODE_BDF(linear_solver=:GMRES),
    space_adapter = RStepAdapter(5, 10, true)
)

adaptiverstepfsp = AdaptiveSparseFsp(
    ode_method = CVODE_BDF(linear_solver=:GMRES),
    space_adapter = SelectiveRStepAdapter(5, 10, true)
)

@btime fspsol1 = solve(model, p0, tspan, fixedrstepfsp, θ, saveat=0.0:20.0:300.0);
@btime fspsol2 = solve(model, p0, tspan, adaptiverstepfsp, θ, saveat=0.0:20.0:300.0);


