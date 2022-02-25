export Propensity, TimeSeparablePropensity, CmeModel 

abstract type Propensity end   

mutable struct TimeSeparablePropensity <: Propensity 
    tfactor::Union{Function,Nothing}
    statefactor::Function         
end

struct CmeModel{PT <: Propensity, IntT <: Integer}
    stoich_matrix::Matrix{IntT}
    propensities::Vector{PT}
end
