# This test suite verifies that DifferentialEquations.jl Exponential integrators can be used with the FSP 
using NumCME
using DifferentialEquations
using StaticArrays 
import Catalyst 

## Test Fixture: telegraph gene expression
Catalyst.@parameters k₀₁ k₁₀ λ γ
bursting_rn = Catalyst.@reaction_network begin 
    k₀₁, G0 --> G1 
    k₁₀, G1 --> G0 
    λ, G1 --> G1 + mRNA 
    γ, mRNA --> ∅
end k₀₁ k₁₀ λ γ

parameter_values = [k₀₁ => 0.05, k₁₀ => 0.1, λ => 5.0, γ => 0.5]
model_from_catalyst = CmeModel(bursting_rn, parameter_values)

fspalgorithm = AdaptiveFspSparse(
    ode_method=LinearExponential(),
    space_adapter=RStepAdapter(5, 10, true)
)
p0 = FspVectorSparse([@MVector [1, 0, 0]], [1.0])
fspsol = solve(model_from_catalyst, p0, (0.0, 100.0), fspalgorithm; saveat=[0.0, 50.0, 100.0])


