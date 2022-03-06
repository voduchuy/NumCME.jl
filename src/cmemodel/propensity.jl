export Propensity, TimeInvariantPropensity, TimeVaryingPropensity, propensity, propensity_timevarying, istimeseparable, istimevarying
export StandardTimeInvariantPropensity, SeparableTimeVaryingPropensity, JointTimeVaryingPropensity
export PropensityGradient, TimeVaryingPropensityGradient, JointTimeVaryingPropensityGradient, SeparableTimeVaryingPropensityGradient, propensitygrad, propensitygrad_timevarying

abstract type Propensity end
abstract type TimeInvariantPropensity <: Propensity end
abstract type TimeVaryingPropensity <: Propensity end

istimevarying(α::Propensity) = false
istimeseparable(α::Propensity) = false

struct StandardTimeInvariantPropensity <: TimeInvariantPropensity
    f::Any
end

struct SeparableTimeVaryingPropensity <: TimeVaryingPropensity
    tfactor::Union{Any,Nothing}
    statefactor::Any
end

struct JointTimeVaryingPropensity <: TimeVaryingPropensity
    f::Any
end

istimevarying(α::TimeVaryingPropensity) = true
istimeseparable(α::SeparableTimeVaryingPropensity) = true

function propensity(f::Any)
    return StandardTimeInvariantPropensity(f)
end

function propensity_timevarying(f::Any)
    return JointTimeVaryingPropensity(f)
end

function propensity_timevarying(xfactor::Any, tfactor::Any)
    return SeparableTimeVaryingPropensity(tfactor, xfactor)
end

# Call methods to make Propensity objects callable
function (a::StandardTimeInvariantPropensity)(x::AbstractVector, p=[])
    a.f(x,p)
end
function (a::StandardTimeInvariantPropensity)(t::AbstractFloat, x::AbstractVector, p=[])
    a.f(x,p)
end
function (a::JointTimeVaryingPropensity)(t::AbstractFloat, x::AbstractVector, p=[])
    a.f(t,x,p)
end
function (a::SeparableTimeVaryingPropensity)(t, x::AbstractVector, p=[])
    a.tfactor(t,p)*a.statefactor(x,p)
end


