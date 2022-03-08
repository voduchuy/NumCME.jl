using Julifsp
using Test 
using BenchmarkTools
using SparseArrays

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
a4 = propensity_timevarying((t,p) -> max(0.0, 1.0-sin(π*t/p[5]))) do x, p
    p[4] * x[3]
end
propensities = [a1, a2, a3, a4]

model = CmeModel(𝕊, propensities, θ)
sensmodel = CmeModelWithSensitivity(model)
space = SparseStateSpace(model.stoich_matrix, [1,0,0])
expand!(space, 200)

sensmat = ForwardSensFspMatrixSparse{Float64}(sensmodel, space)

pcount = get_parameter_count(model)
n = get_rowcount(sensmat.fspmatrix)
v = ones(n + pcount*n)
v[1] = 1.0
out = similar(v)
matvec!(out, 10.0, sensmat, v)
@test sum(out) == 0.0