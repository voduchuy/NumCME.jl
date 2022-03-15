using NumCME
using Test 
using BenchmarkTools
using SparseArrays 
using ForwardDiff

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
a4 = propensity((t,p) -> max(0.0, 1.0-sin(π*t/p[5]))) do x, p
    p[4] * x[3]
end
propensities = [a1, a2, a3, a4]
sensmodel = CmeModelWithSensitivity(CmeModel(𝕊, propensities, θ))

# Test correctness of the generated propensity gradients across multiple CME states and times
test_space = StateSpaceSparse(𝕊, x₀)
expand!(test_space, 100)
test_states = get_states(test_space)
test_times = 0.0:1.0:40.0

a1grad(x, p) = [x[1], 0.0, 0.0, 0.0, 0.0]
a2grad(x, p) = [0.0, x[2], 0.0, 0.0, 0.0]
a3grad(x, p) = [0.0, 0.0, x[2], 0.0, 0.0]
a4grad(t, x, p) = (1.0-sin(π*t/p[5]) > 0) ? [0.0, 0.0, 0.0, (1.0-sin(π*t/p[5]))*x[3], p[4]*x[3]*cos(π*t/p[5])*π*t/(p[5]^2)] : [0.0, 0.0, 0.0, 0.0, 0.0]
propensity_grads = get_propensity_gradients(sensmodel)

@test prod([a1grad(x, θ) ≈ propensity_grads[1](x, θ) for x in test_states])
@test prod([a2grad(x, θ) ≈ propensity_grads[2](x, θ) for x in test_states])
@test prod([a3grad(x, θ) ≈ propensity_grads[3](x, θ) for x in test_states])

tmp = []
for t in test_times 
    push!(tmp, prod([isapprox(a4grad(t, x, θ),propensity_grads[4](t, x, θ), rtol=1.0E-8, atol=10*eps()) for x in test_states]))
    for x in test_states 
        !(isapprox(a4grad(t, x, θ),propensity_grads[4](t, x, θ), rtol=1.0E-8, atol=10*eps())) && println(a4grad(t, x, θ)," ",propensity_grads[4](t, x, θ))
    end
end
@test prod(tmp)

# Test correctness of auto-determined sparsity pattern of propensity gradients 
correct_sparsity = sparse([1, 2, 3, 4, 4],[1, 2, 3, 4, 5],[true, true, true, true, true])
sparsity_patterns = get_gradient_sparsity_patterns(sensmodel)
@test correct_sparsity == sparsity_patterns












