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

θ = [k₀₁, k₁₀, λ, γ]

α₁ = propensity() do x, p 
    p[1]*x[1]
end
α₂ = propensity() do x, p 
    p[2]*x[2]
end
α₂tv = propensity((t,p) -> max(0.0, 1.0 - sin(π * t / 2))) do x, p 
    p[2]*x[2]
end
α₂tvj = propensity() do t, x, p 
    max(0.0, 1.0 - sin(π * t / 2))*p[2]*x[2]
end
α₃ = propensity() do x, p 
    p[3]*x[2]
end
α₄ = propensity() do x, p 
    p[4]*x[3]
end

propensities_ti = [α₁, α₂, α₃, α₄]
propensities_tv = [α₁, α₂tv, α₃, α₄]
propensities_tvj = [α₁, α₂tvj, α₃, α₄]

𝔛 = SparseStateSpace(𝕊, x₀)
expand!(𝔛, 2)
𝐀 = FspMatrixSparse(𝔛, propensities_ti, parameters=θ)
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
A1 = FspMatrixSparse(𝔛, propensities_tv, parameters=θ)
@test size(A1, 1) == get_state_count(𝔛) + get_sink_count(𝔛)
@test size(A1, 2) == get_state_count(𝔛) + get_sink_count(𝔛)
𝐯 = ones(Float64, size(A1, 1))
w1 = A1(1.0) * 𝐯
@test sum(w1) ≈ 0.0 atol = 1.0e-14
w1 = A1 * 𝐯
@test sum(w1) ≈ 0.0 atol = 1.0e-14

A2 = FspMatrixSparse(𝔛, propensities_tvj, parameters=θ)
w2 = A2(1.0) * 𝐯
@test sum(w2) ≈ 0.0 atol = 1.0e-14
w2 = A2 * 𝐯

@test sum(w2) ≈ 0.0 atol = 1.0e-14
@test norm(w1 -w2) ≈ 0


