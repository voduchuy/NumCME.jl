export StandardTimeInvariantPropensityForwardDiff, JointTimeVaryingPropensityForwardDiff, SeparableTimeVaryingPropensityForwardDiff, tfactorgrad, statefactorgrad

struct StandardTimeInvariantPropensityForwardDiff{TIP <: TimeInvariantPropensity} <: TimeInvariantPropensityGradient
    propensity::TIP 
    cfg::ForwardDiff.GradientConfig    
    species_count::Integer 
    parameter_count::Integer
end
function StandardTimeInvariantPropensityForwardDiff(f::TIP, species_count, parameter_count) where {TIP <: TimeInvariantPropensity}
    return StandardTimeInvariantPropensityForwardDiff(f, ForwardDiff.GradientConfig(nothing, zeros(species_count+parameter_count)), species_count, parameter_count)
end
function (g::StandardTimeInvariantPropensityForwardDiff)(x::AbstractVector, p::AbstractVector)
    return ForwardDiff.gradient(z->g.propensity(z[1:g.species_count],z[g.species_count+1:end]), [x;p], g.cfg)[end-g.parameter_count+1:end]
end

struct JointTimeVaryingPropensityForwardDiff{TVP <: JointTimeVaryingPropensity} <: TimeVaryingPropensityGradient
    species_count::Integer 
    parameter_count::Integer
    propensity::TVP 
    cfg::ForwardDiff.GradientConfig    
end
function JointTimeVaryingPropensityForwardDiff(f::TVP, species_count, parameter_count) where {TVP <: JointTimeVaryingPropensity}
    return JointTimeVaryingPropensityForwardDiff(f, ForwardDiff.GradientConfig(nothing, zeros(species_count+parameter_count)), species_count, parameter_count)
end
function (g::JointTimeVaryingPropensityForwardDiff)(t::AbstractFloat, x::Vector, p::Vector)
    return ForwardDiff.gradient(z->g.propensity(z[1], z[1:g.species_count],z[g.species_count+1:end]), [t;x;p], g.cfg)[end-g.parameter_count+1:end]
end

struct SeparableTimeVaryingPropensityForwardDiff{TVP <:SeparableTimeVaryingPropensity} <: TimeVaryingPropensityGradient 
    species_count::Integer 
    parameter_count::Integer
    propensity::TVP 
    tfactorcfg::ForwardDiff.GradientConfig
    statefactorcfg::ForwardDiff.GradientConfig    
end
function SeparableTimeVaryingPropensityForwardDiff(f::TVP, species_count, parameter_count) where {TVP <: SeparableTimeVaryingPropensity}
    return SeparableTimeVaryingPropensityForwardDiff(species_count, parameter_count, f, ForwardDiff.GradientConfig(nothing, zeros(1+parameter_count)), ForwardDiff.GradientConfig(nothing, zeros(species_count+parameter_count)))
end
function tfactorgrad(g::SeparableTimeVaryingPropensityForwardDiff, t::AbstractFloat, p::AbstractVector)
    return ForwardDiff.gradient(
        z->g.propensity.tfactor(z[1], z[2:end]), [t;p], g.tfactorcfg
    )[end-g.parameter_count+1:end]
end
function statefactorgrad(g::SeparableTimeVaryingPropensityForwardDiff, x::AbstractVector, p::AbstractVector)
    return ForwardDiff.gradient(
        z->g.propensity.statefactor(z[1:g.species_count], z[g.species_count+1:end]), [x;p], g.statefactorcfg
    )[end-g.parameter_count+1:end]
end
function (g::SeparableTimeVaryingPropensityForwardDiff)(t::AbstractFloat, x::AbstractVector, p::AbstractVector)
    return g.propensity.tfactor(t, p)*statefactorgrad(g, x, p) + g.propensity.statefactor(x, p)*tfactorgrad(g, t, p)
end

propensity_forwarddiff(f::SeparableTimeVaryingPropensity, species_count, parameter_count) = SeparableTimeVaryingPropensityForwardDiff(f, species_count, parameter_count)
propensity_forwarddiff(f::JointTimeVaryingPropensity, species_count, parameter_count) = JointTimeVaryingPropensityForwardDiff(f, species_count, parameter_count)
propensity_forwarddiff(f::StandardTimeInvariantPropensity, species_count, parameter_count) = StandardTimeInvariantPropensityForwardDiff(f, species_count, parameter_count)