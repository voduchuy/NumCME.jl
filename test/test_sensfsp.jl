using NumCME
using Test
using Sundials: CVODE_BDF

𝕊 = [[-1, 1, 0] [1, -1, 0] [0, 0, 1] [0, 0, -1]]
x₀ = [1, 0, 0]
k₀₁ = 0.05
k₁₀ = 0.1
λ = 5.0
γ = 0.5
L = 20.0
θ = [k₀₁, k₁₀, λ, γ, L]
a1 = propensity() do x, p
    p[1] * x[1]
end
a2 = propensity() do x, p
    p[2] * x[2]
end
a3 = propensity() do x, p
    p[3] * x[2]
end
a4 = propensity((t, p) -> max(0.0, 1.0 - sin(π * t / p[5]))) do x, p
    p[4] * x[3]
end
propensities = [a1, a2, a3, a4]
sensmodel = CmeModelWithSensitivity(CmeModel(𝕊, propensities, θ))
init_cond = forwardsens_initial_condition([x₀], [1.0], [[0.0] for i in 1:5])
sensfspalgorithm = AdaptiveForwardSensFspSparse(
    ode_method = CVODE_BDF(linear_solver = :GMRES),
    space_adapter = ForwardSensRStepAdapter(10, 10, true)
)
@time senssol = solve(sensmodel,
    init_cond,
    (0.0, 400.0),
    sensfspalgorithm;
    saveat = [],
    fsptol = 1.0E-6,
    odeatol = 1.0E-10,
    odertol = 1.0E-4)





