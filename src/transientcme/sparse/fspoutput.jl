export FspOutputSparse, FspOutputSliceSparse

"""
    $(TYPEDEF)

Struct to store Finite State Projection outputs based on the sparse representation of the FSP solution.

# Fields 
- `t`: Array of solution output times.
- `p`: Array of `FspVectorSparse` instances. `p[i]` is the solution at time `t[i]`.
- `sinks`: Probability mass accumulated at the sink states over time. 

# Usage 
If `sol` is of type `FspOutputSparse`, `sol[i]` will return a slice, of type `FspOutputSliceSparse` of the solution set at the `i`-th index. 

# See also 
[`FspOutputSliceSparse`](@ref), [`length(::FspOutputSparse)`](@ref).
"""
Base.@kwdef struct FspOutputSparse{NS,IntT<:Integer,RealT<:AbstractFloat}
    t::Vector{RealT}
    p::Vector{FspVectorSparse{NS,IntT,RealT}}
    sinks::Vector{Vector{RealT}}
end

"""
    $(TYPEDEF)

Struct to store the FSP solution at a single time. Fields: `t`, `p`, `sinks`.
"""
struct FspOutputSliceSparse
    t::AbstractFloat 
    p::FspVectorSparse
    sinks::Vector
end

function FspOutputSparse{NS}() where {NS, IntT <: Integer, RealT <: AbstractFloat}
    return FspOutputSparse{NS, Int64, Float64}(
        t = Vector{Float64}(),
        p = Vector{FspVectorSparse{NS, IntT, RealT}}(),
        sinks = Vector{Vector{RealT}}()
    )
end

function Base.getindex(fspoutput::FspOutputSparse, ind::Integer)
    if ind > length(fspoutput.t)
        throw(ArgumentError("Requested index exceeds array limit."))
    end
    return FspOutputSliceSparse(
        fspoutput.t[ind],
        fspoutput.p[ind],
        fspoutput.sinks[ind]
    )
    return 
end

function Base.getindex(fspoutput::FspOutputSparse, inds::Vector{T}) where {T <: Integer}    
    return [fspoutput[i] for i in inds]    
end

function Base.lastindex(fspoutput::FspOutputSparse)
    return lastindex(fspoutput.t)
end

"""
$(TYPEDSIGNATURES)

Return number of FSP solutions in the solution set.
"""
function Base.length(fspoutput::FspOutputSparse)
    return length(fspoutput.t)
end
function Base.size(fspoutput::FspOutputSparse)
    return size(fspoutput.t)
end

using Printf 
function Base.show(io::IO, outputslice::FspOutputSliceSparse)
    println(io, "FSP solution output slice:")
    @printf(io, "t = %.2e \n", outputslice.t)
    println(io, "sinks = $(outputslice.sinks)")
    @printf(io, "p = \n")
    show(outputslice.p)    
end