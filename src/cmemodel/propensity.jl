export Propensity, TimeInvariantPropensity, TimeVaryingPropensity, propensity, propensity, istimeseparable, istimevarying
export StandardTimeInvariantPropensity, SeparableTimeVaryingPropensity, JointTimeVaryingPropensity
export PropensityGradient, TimeVaryingPropensityGradient, JointTimeVaryingPropensityGradient, SeparableTimeVaryingPropensityGradient, propensitygrad, propensitygrad_timevarying

"""
Base type for all propensities.
"""
abstract type Propensity end

abstract type TimeInvariantPropensity <: Propensity end
abstract type TimeVaryingPropensity <: Propensity end

istimevarying(α::Propensity) = false
istimeseparable(α::Propensity) = false

"""
    StandardTimeInvariantPropensity

Propensity function for a reaction with rate not dependent on time. 

# Examples

A standard time-invariant propensity can be created by wrapping around a Julia function. The function must accept two inputs `x`, `p` where `x` is the CME state vector and `p` is the vector of parameters.
```@repl
julia> f(x, p) = p[1]*x[1]
f (generic function with 1 method)
julia> α = StandardTimeInvariantPropensity(f)
StandardTimeInvariantPropensity(f)
julia> α([1,2], [0.1,0.2])
0.1
```

One can also input `f` into the `propensity()` function. Provided that `f` is a function of two arguments with the CME state argument followed by the parameter argument, `propensity()` will output a standard time-invariant propensity. Continuing the example above, 
```@repl
julia> β = propensity(f)
StandardTimeInvariantPropensity(f)
julia> β([1,2], [0.1,0.2]) == α([1,2], [0.1,0.2])
true 
```

One can also take advantage of Julia's `do` block to make the code more concise.
```@repl
julia> γ = propensity() do x,p
        x[1]*p[1]
    end 
julia> γ([1,2], [0.1,0.2])
0.1
```
"""
struct StandardTimeInvariantPropensity <: TimeInvariantPropensity
    f::Any
end

"""
    SeparableTimeVaryingPropensity

Reaction propensity with a time-varying rate that can be factorized into a function that depends only on the time variable (and parameters) and another function that does not depend on time. Subsequent types that build upon this type (such as `FspMatrixSparse`) takes advantage of time-space separability to avoid repeating expensive computations.

# Examples

A separable time-varying propensity can be created by specifying two Julia functions. The first function represents the time-varying factor and accepts two input arguments `t` (for time) and `p` (for parameters), the second function represents the time-invariant factor and accepts two input arguments `x` (for CME state) and `p` (for parameters).
```@repl
julia> c(t,p) = (1.0+cos(π*t/p[2]))
c (generic function with 1 method)
julia> f(x, p) = p[1]*x[1]
f (generic function with 1 method)
julia> α = SeparableTimeVaryingPropensity(c,f)
SeparableTimeVaryingPropensity(c,f)
julia> α(20.0,[1,2], [0.1,0.2])
0.2
```

One can also use a method of the `propensity()` function to create a `SeparableTimeVaryingPropensity` instance. Provided that `f` is a function of two arguments with the CME state argument followed by the parameter argument, `propensity(f,c)` (note that the time-varying factor comes _after_ the time-invariant factor) will output a standard time-invariant propensity. Continuing the example above, 
```@repl
julia> β = propensity(f,c)
SeparableTimeVaryingPropensity(f, c)
julia> β(20.0, [1,2], [0.1,0.2]) == α(20.0, [1,2], [0.1,0.2])
true 
```

One can also take advantage of Julia's `do` block to make the code more concise.
```@repl
julia> γ = propensity((1.0+cos(π*t/p[2]))) do x,p
        x[1]*p[1]
    end 
julia> SeparableTimeVaryingPropensity(var"#19#21"(), var"#18#20"())
julia> γ(20.0, [1,2], [0.1,0.2])
0.2
```
"""
struct SeparableTimeVaryingPropensity <: TimeVaryingPropensity
    tfactor::Union{Any,Nothing}
    statefactor::Any
end

"""
    JointTimeVaryingPropensity

Reaction propensity with a time-varying rate. For speed, we recommend using `SeparableTimeVaryingPropensity` if the propensity function has a separable structure.

# Examples

A time-varying propensity can be created by wrapping a function with three input arguments, `t` (for time) `x` (for CME state) and `p` (for parameters).

```@repl
julia> f(t, x, p) = (1.0+cos(π*t/p[2]))*p[1]*x[1]
f (generic function with 1 method)
julia> α = JointTimeVaryingPropensity(c,f)
JointTimeVaryingPropensity(c,f)
julia> α(20.0,[1,2], [0.1,0.2])
0.2
```

One can also use a method of the `propensity()` function to create a `JointTimeVaryingPropensity` instance. Provided that `f` is a function of three arguments with time as the first argument, the CME state as the second argument followed by the parameter argument.

```@repl
julia> β = propensity(f)
julia> β(20.0, [1,2], [0.1,0.2]) == α(20.0, [1,2], [0.1,0.2])
true 
```
"""
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


