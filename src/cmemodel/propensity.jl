export Propensity, TimeInvariantPropensity, TimeVaryingPropensity, propensity, propensity, istimeseparable, istimevarying
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
    m = methods(f)[1]
    if m.nargs == 3
        return StandardTimeInvariantPropensity(f)
    elseif m.nargs == 4
        return JointTimeVaryingPropensity(f)
    else 
        throw(ArgumentError("The callable passed to `propensity()` must have either two arguments (x,p) or three arguments (t,x,p)."))
    end
end

function propensity(xfactor::Any, tfactor::Any)
    return SeparableTimeVaryingPropensity(tfactor, xfactor)
end

# Call methods to make Propensity objects callable
function (a::StandardTimeInvariantPropensity)(x::AbstractVector, p=[])
    a.f(x,p)
end
function (a::JointTimeVaryingPropensity)(t::Real, x::AbstractVector, p=[])
    a.f(t,x,p)
end
function (a::SeparableTimeVaryingPropensity)(t::Real, x::AbstractVector, p=[])
    a.tfactor(t,p)*a.statefactor(x,p)
end


