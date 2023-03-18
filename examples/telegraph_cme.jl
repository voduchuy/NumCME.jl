using NumCME
using BenchmarkTools
using Catalyst
using Sundials: CVODE_BDF
using StaticArrays: @MVector


fspalgorithm = AdaptiveFspSparse(
    ode_method=CVODE_BDF(linear_solver=:GMRES),
    space_adapter=RStepAdapter(5, 10, true)
)


# Bursting gene model definition using basic Julia
𝕊 = [[-1, 1, 0] [1, -1, 0] [0, 0, 1] [0, 0, -1]]
x₀ = [1, 0, 0]

a1 = propensity() do x, p
    p[1] * x[1]
end
a2 = propensity() do x, p
    p[2] * x[2]
end
a3 = propensity() do x, p
    p[3] * x[2]
end
a4 = propensity() do x, p
    p[4] * x[3]
end

k₀₁ = 0.05
k₁₀ = 0.1
λ = 5.0
γ = 0.5
θ = [k₀₁, k₁₀, λ, γ]

model = CmeModel(𝕊, [a1, a2, a3, a4], θ)

p0 = FspVectorSparse([@MVector [1, 0, 0]], [1.0])
tspan = (0.0, 300.0)

fspsol1 = solve(model, p0, tspan, fspalgorithm)

# Bursting model definition using Catalyst 
bursting_rn = @reaction_network begin 
    k₀₁, G0 --> G1 
    k₁₀, G1 --> G0 
    λ, G1 --> G1 + mRNA 
    γ, mRNA --> ∅
end 

parameter_values = [k₀₁ => 0.05, k₁₀ => 0.1, λ => 5.0, γ => 0.5]
model_from_catalyst = CmeModel(bursting_rn, parameter_values)

fspsol2 = solve(model_from_catalyst, p0, tspan, fspalgorithm)

# Check that the two ways to code the model lead to the same numerical outputs 
@assert length(fspsol1) == length(fspsol2)
for i ∈ 1:length(fspsol1)
    @assert get_states(fspsol1[i].p) == get_states(fspsol2[i].p)
    @assert get_values(fspsol1[i].p) == get_values(fspsol2[i].p)
end

# Which method is faster?
@btime solve(model, p0, tspan, fspalgorithm);
@btime solve(model_from_catalyst, p0, tspan, fspalgorithm);





