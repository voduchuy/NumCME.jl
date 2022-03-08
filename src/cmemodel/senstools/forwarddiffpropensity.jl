export propensity_forwarddiff

function propensity_forwarddiff(f::StandardTimeInvariantPropensity, parameter_count)
    pardiffs = []
    for ip in 1:parameter_count
        push!(pardiffs, 
            (x, p) -> ForwardDiff.derivative(z -> f(x,[p[1:ip-1];z;p[ip+1:end]]), p[ip])
        )
    end
    return StandardTimeInvariantPropensityGradient(pardiffs)
end

function propensity_forwarddiff(f::JointTimeVaryingPropensity, parameter_count)
    pardiffs = []
    for ip in 1:parameter_count
        push!(pardiffs, 
            (t, x, p) -> ForwardDiff.derivative(z -> f(t, x,[p[1:ip-1];z;p[ip+1:end]]), p[ip])
        )
    end
    return JointTimeVaryingPropensityGradient(pardiffs)
end

function propensity_forwarddiff(f::SeparableTimeVaryingPropensity, parameter_count)
    tfactor_pardiffs = []
    for ip in 1:parameter_count
        push!(tfactor_pardiffs, 
            (t, p) -> ForwardDiff.derivative(z -> f.tfactor(t,[p[1:ip-1];z;p[ip+1:end]]), p[ip])
        )
    end
    statefactor_pardiffs = []
    for ip in 1:parameter_count
        push!(statefactor_pardiffs, 
            (x, p) -> ForwardDiff.derivative(z -> f.statefactor(x,[p[1:ip-1];z;p[ip+1:end]]), p[ip])
        )
    end
    return SeparableTimeVaryingPropensityGradient(f.tfactor, f.statefactor, tfactor_pardiffs, statefactor_pardiffs)
end