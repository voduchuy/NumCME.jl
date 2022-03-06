export PropensityGradient, TimeInvariantPropensityGradient, TimeVaryingPropensityGradient, statefactorgrad, tfactorgrad

abstract type PropensityGradient end
abstract type TimeInvariantPropensityGradient <: PropensityGradient end
abstract type TimeVaryingPropensityGradient <: PropensityGradient end

struct StandardTimeInvariantPropensityGradient <: TimeInvariantPropensityGradient
    gradf::Any
end

struct SeparableTimeVaryingPropensityGradient <: TimeVaryingPropensityGradient
    tfactor::Any
    statefactor::Any
    gradtfactor::Union{Any,Nothing}
    gradstatefactor::Any
end

struct JointTimeVaryingPropensityGradient <: TimeVaryingPropensityGradient
    gradf::Any
end

function propensitygrad(∇f::Any)
    return StandardTimeInvariantPropensityGradient(∇f)
end

function propensitygrad_timevarying(∇f::Any)
    return JointTimeVaryingPropensityGradient(∇f)
end

function propensitygrad_timevarying(tfactor::Any, statefactor::Any, ∇tfactor::Any, ∇statefactor::Any)
    return SeparableTimeVaryingPropensityGradient(tfactor, statefactor, ∇tfactor, ∇statefactor)
end

(g::StandardTimeInvariantPropensityGradient)(x::AbstractVector, p::AbstractVector) = g.gradf(x, p)
(g::JointTimeVaryingPropensityGradient)(t::AbstractFloat, x::AbstractVector, p::AbstractVector) = g.gradf(t, x, p)
(g::SeparableTimeVaryingPropensityGradient)(t::AbstractFloat, x::AbstractVector, p::AbstractVector) = g.tfactor(t,p)*g.gradtfactor(t,p) + g.statefactor(x,p)*g.gradstatefactor(x,p)

include("senstools/sparsity_pattern.jl")
include("senstools/forwarddiffpropensity.jl")