export PropensityGradient, TimeInvariantPropensityGradient, TimeVaryingPropensityGradient, statefactorgrad, tfactorgrad

"""
Base type for all propensity gradients.
"""
abstract type PropensityGradient end
abstract type TimeInvariantPropensityGradient <: PropensityGradient end
abstract type TimeVaryingPropensityGradient <: PropensityGradient end

"""
Store the partial derivatives of a `StandardTimeInvariantPropensity` instance. 

# Fields 
- `pardiffs`: Vector of callables, whose length is the number of model parameters. Each of these callables has syntax `pardiffs[i](x,p)::Real` where `x` is the state vector and `p` the parameter vector.
"""
struct StandardTimeInvariantPropensityGradient <: TimeInvariantPropensityGradient
    pardiffs::Vector{Any}
end

"""
Store the partial derivatives of a `SeparableTimeVaryingPropensity` instance. 

# Fields 
- `tfactor`: Callable in the form `tfactor(t,p)` where `t` is time, `p` is parameter vector.
- `statefactor`: Callable in the form `statefactor(x,p)` where `x` is CME state, `p` the parameter vector.
- `tfactor_pardiffs`: Vector of callables. Each of these callables has syntax `tfactor_pardiffs[i](t,p)::Real` where `t` is time variable and `p` the parameter vector.
- `statefactor_pardiffs`: Vector of callables, whose length is the number of model parameters. Each of these callables has syntax `statefactor_pardiffs[i](x,p)::Real` where `x` is the state vector and `p` the parameter vector.
"""
struct SeparableTimeVaryingPropensityGradient <: TimeVaryingPropensityGradient
    tfactor::Any
    statefactor::Any
    tfactor_pardiffs::Vector{Any}
    statefactor_pardiffs::Vector{Any}
end

"""
Store the partial derivatives of a `JointTimeVaryingPropensity` instance. 

# Fields 
- `pardiffs`: Vector of callables, whose length is the number of model parameters. Each of these callables has syntax `pardiffs[i](t,x,p)::Real` where `t` is time, `x` is the state vector and `p` the parameter vector.
"""
struct JointTimeVaryingPropensityGradient <: TimeVaryingPropensityGradient
    pardiffs::Vector{Any}
end

function propensitygrad(∇f::Vector{Any})
    return StandardTimeInvariantPropensityGradient(∇f)
end
function propensitygrad_timevarying(∇f::Vector{Any})
    return JointTimeVaryingPropensityGradient(∇f)
end
function propensitygrad_timevarying(tfactor::Any, statefactor::Any, ∇tfactor::Vector{Any}, ∇statefactor::Vector{Any})
    return SeparableTimeVaryingPropensityGradient(tfactor, statefactor, ∇tfactor, ∇statefactor)
end

function (g::StandardTimeInvariantPropensityGradient)(x::AbstractVector, p::AbstractVector)
    return [g.pardiffs[i](x, p) for i in 1:length(g.pardiffs)]
end
function (g::JointTimeVaryingPropensityGradient)(t::AbstractFloat, x::AbstractVector, p::AbstractVector)
    return [g.pardiffs[i](t, x, p) for i in 1:length(g.pardiffs)]
end
function (g::SeparableTimeVaryingPropensityGradient)(t::AbstractFloat, x::AbstractVector, p::AbstractVector)
    return [
        g.statefactor(x, p) * g.tfactor_pardiffs[i](t, p) + g.tfactor(t, p) * g.statefactor_pardiffs[i](x, p)
        for i in 1:length(g.tfactor_pardiffs)
    ]
end

get_pardiffs(g::Union{StandardTimeInvariantPropensityGradient,JointTimeVaryingPropensityGradient}) = g.pardiffs 
get_single_pardiff(g::Union{StandardTimeInvariantPropensityGradient, JointTimeVaryingPropensityGradient}, ip) = g.pardiffs[ip]
get_tfactor_pardiffs(g::SeparableTimeVaryingPropensityGradient) = g.tfactor_pardiffs
get_statefactor_pardiffs(g::SeparableTimeVaryingPropensityGradient) = g.statefactor_pardiffs
get_single_tfactor_pardiff(g::SeparableTimeVaryingPropensityGradient, ip) = g.tfactor_pardiffs[ip]
get_single_statefactor_pardiff(g::SeparableTimeVaryingPropensityGradient, ip) = g.statefactor_pardiffs[ip]


include("senstools/sparsity_pattern.jl")
include("senstools/forwarddiffpropensity.jl")