using Julifsp
using Test 
using BenchmarkTools
using SparseArrays

𝕊 = [[1] [-1]]
λ = 10.0
γ = 5.0
a1 = propensity() do x, p 
    p[1]
end
a2 = propensity() do x, p 
    p[2]
end
sensmodel = CmeModelWithSensitivity(CmeModel(𝕊, [a1, a2], [λ, γ]))
test_space = StateSpaceSparse(𝕊, [[i] for i in 1:100])
sensmat = ForwardSensFspMatrixSparse{Float64}(sensmodel, test_space)
v = ones(3*get_rowcount(sensmat.fspmatrix))
out = similar(v)
matvec!(out, 0.0, sensmat, v)
@test sum(out) ≈ 0.0