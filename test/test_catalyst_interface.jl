using NumCME
using Catalyst
using StaticArrays: MVector, @MVector 
using Test

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

cmemodel1 = CmeModel(𝕊, [a1,a2,a3,a4], θ)

rn = @reaction_network begin
    k01, G0 --> G1
    k10, G1 --> G0
    α, G1 --> G1 + RNA
    γ*max(0.0, 1.0 - sin(π * t / L)), RNA --> ∅
end k01 k10 α γ L

cmemodel2 = CmeModel(rn, θ)
@test get_species_count(cmemodel2) == 3
@test get_reaction_count(cmemodel2) == 4
@test get_parameter_count(cmemodel2) == 5
for i in 1:3
    @test !istimevarying(cmemodel2.propensities[i])
end
@test istimevarying(cmemodel2.propensities[4])

test_times = [0.0, 10.0, 20.0, 30.0]
test_space = StateSpaceSparse(𝕊, [1,0,0])
expand!(test_space, 100)
pass = true 
for t in test_times 
    for state in get_states(test_space)
        global pass &= cmemodel1.propensities[1](state,θ) ≈ cmemodel2.propensities[1](state,θ)
        pass &= cmemodel1.propensities[2](state,θ) ≈ cmemodel2.propensities[2](state,θ)
        pass &= cmemodel1.propensities[3](state,θ) ≈ cmemodel2.propensities[3](state,θ)
        pass &= cmemodel1.propensities[4](t,state,θ) ≈ cmemodel2.propensities[4](t,state,θ)
        (!pass) && break
    end
end
@test pass 







