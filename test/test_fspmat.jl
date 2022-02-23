using EtFsp.statespace: FspStateSpaceBasic, rstep_expand!, FspStateSpace, get_state_count, get_sink_count 
using EtFsp.fspmatrix: Propensity, FspMatrix
using Test

## Bursting gene model 
𝕊 = [[-1, 1, 0] [1, -1, 0] [0, 0, 1] [0, 0, -1]]
x₀ = [1,0,0]
k₀₁ = 0.05
k₁₀ = 0.1
λ = 5.0
γ = 1.0
𝕻 = [
    Propensity(nothing, (G₀::Int, G₁::Int, RNA::Int)->k₀₁*G₀)
    Propensity(nothing, (G₀::Int, G₁::Int, RNA::Int)->k₁₀*G₁)
    Propensity(nothing, (G₀::Int, G₁::Int, RNA::Int)->λ*G₁)
    Propensity(nothing, (G₀::Int, G₁::Int, RNA::Int)->γ*RNA)
]
𝔛 = FspStateSpaceBasic(𝕊, x₀)
rstep_expand!(𝔛, 2)
𝐀 = FspMatrix(𝔛, 𝕻)
@test size(𝐀, 1) == get_state_count(𝔛) + get_sink_count(𝔛)
@test size(𝐀, 2) == get_state_count(𝔛) + get_sink_count(𝔛)
𝐯 = ones(Float64, size(𝐀, 1))
𝐰 = 𝐀(1.0)*𝐯 
@test sum(𝐰) ≈ 0.0 atol=1.0e-14 
𝐰 = 𝐀*𝐯
@test sum(𝐰) ≈ 0.0 atol=1.0e-14 

## Bursting gene model with time-varying rate 
𝕊 = [[-1, 1, 0] [1, -1, 0] [0, 0, 1] [0, 0, -1]]
x₀ = [1,0,0]
k₀₁ = 0.05
k₁₀ = 0.1
λ = 5.0
γ = 1.0
𝕻 = [
    Propensity(nothing, (G₀::Int, G₁::Int, RNA::Int)->k₀₁*G₀)
    Propensity(t->max(0.0, 1.0 - sin(π*t/2)), (G₀::Int, G₁::Int, RNA::Int)->k₁₀*G₁)
    Propensity(nothing, (G₀::Int, G₁::Int, RNA::Int)->λ*G₁)
    Propensity(nothing, (G₀::Int, G₁::Int, RNA::Int)->γ*RNA)
]
𝔛 = FspStateSpaceBasic(𝕊, x₀)
rstep_expand!(𝔛, 2)
𝐀 = FspMatrix(𝔛, 𝕻)
@test size(𝐀, 1) == get_state_count(𝔛) + get_sink_count(𝔛)
@test size(𝐀, 2) == get_state_count(𝔛) + get_sink_count(𝔛)
𝐯 = ones(Float64, size(𝐀, 1))
𝐰 = 𝐀(1.0)*𝐯 
@test sum(𝐰) ≈ 0.0 atol=1.0e-14 
𝐰 = 𝐀*𝐯 
@test sum(𝐰) ≈ 0.0 atol=1.0e-14 

## Test the integration of the FSP system 
𝕊 = [[-1, 1, 0] [1, -1, 0] [0, 0, 1] [0, 0, -1]]
x₀ = [1,0,0]
k₀₁ = 0.05
k₁₀ = 0.1
λ = 5.0
γ = 1.0
𝕻 = [
    Propensity(nothing, (G₀::Int, G₁::Int, RNA::Int)->k₀₁*G₀)
    Propensity(t->max(0.0, 1.0 - sin(π*t/2)), (G₀::Int, G₁::Int, RNA::Int)->k₁₀*G₁)
    Propensity(nothing, (G₀::Int, G₁::Int, RNA::Int)->λ*G₁)
    Propensity(nothing, (G₀::Int, G₁::Int, RNA::Int)->γ*RNA)
]
𝔛 = FspStateSpaceBasic(𝕊, x₀)
rstep_expand!(𝔛, 10)
𝐀 = FspMatrix(𝔛, 𝕻)
𝐩₀ = zeros(Float64, size(𝐀, 1))
