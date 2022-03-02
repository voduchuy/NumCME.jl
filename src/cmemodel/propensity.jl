export Propensity, TimeInvariantPropensity, TimeVaryingPropensity, propensity, propensity_timevarying, istimeseparable, istimevarying
export StandardTimeIndependentPropensity, SeparableTimeVaryingPropensity, JointTimeVaryingPropensity
export PropensityGradient, TimeVaryingPropensityGradient, JointTimeVaryingPropensityGradient, SeparableTimeVaryingPropensityGradient, propensitygrad, propensitygrad_timevarying

abstract type Propensity end
abstract type TimeInvariantPropensity <: Propensity end
abstract type TimeVaryingPropensity <: Propensity end

istimevarying(α::Propensity) = false
istimeseparable(α::Propensity) = false

struct StandardTimeIndependentPropensity <: TimeInvariantPropensity
    f::Function
end

struct SeparableTimeVaryingPropensity <: TimeVaryingPropensity
    tfactor::Union{Function,Nothing}
    statefactor::Function
end

struct JointTimeVaryingPropensity <: TimeVaryingPropensity
    f::Function
end

istimevarying(α::TimeVaryingPropensity) = true
istimeseparable(α::SeparableTimeVaryingPropensity) = true

function propensity(f::Function)
    return StandardTimeIndependentPropensity(f)
end

function propensity_timevarying(f::Function)
    return JointTimeVaryingPropensity(f)
end

function propensity_timevarying(tfactor::Function, xfactor::Function)
    return SeparableTimeVaryingPropensity(tfactor, xfactor)
end

abstract type PropensityGradient end
abstract type TimeInvariantPropensityGradient <: PropensityGradient end
abstract type TimeVaryingPropensityGradient <: PropensityGradient end

struct StandardTimeIndependentPropensityGradient <: TimeInvariantPropensityGradient
    gradf::Function
end

struct SeparableTimeVaryingPropensityGradient <: TimeVaryingPropensityGradient
    tfactor::Function
    statefactor::Function
    gradtfactor::Union{Function,Nothing}
    gradstatefactor::Function
end

struct JointTimeVaryingPropensityGradient <: TimeVaryingPropensityGradient
    gradf::Function
end

function propensitygrad(∇f::Function)
    return StandardTimeIndependentPropensityGradient(∇f)
end

function propensitygrad_timevarying(∇f::Function)
    return JointTimeVaryingPropensityGradient(∇f)
end

function propensitygrad_timevarying(tfactor::Function, statefactor::Function, ∇tfactor::Function, ∇xfactor::Function)
    return SeparableTimeVaryingPropensityGradient(tfactor, statefactor, ∇tfactor, ∇xfactor)
end
