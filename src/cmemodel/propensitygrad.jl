export PropensityGradient, TimeInvariantPropensityGradient, TimeVaryingPropensityGradient, statefactorgrad, tfactorgrad

abstract type PropensityGradient end
abstract type TimeInvariantPropensityGradient <: PropensityGradient end
abstract type TimeVaryingPropensityGradient <: PropensityGradient end

struct StandardTimeInvariantPropensityGradient <: TimeInvariantPropensityGradient
    pardiffs::Vector{Any}
end


struct SeparableTimeVaryingPropensityGradient <: TimeVaryingPropensityGradient
    tfactor::Any
    statefactor::Any
    tfactor_pardiffs::Vector{Any}
    statefactor_pardiffs::Vector{Any}
end

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