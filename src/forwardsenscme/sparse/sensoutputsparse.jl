"""
    ForwardSensFspOutputSparse{NS,IntT<:Integer,RealT<:AbstractFloat}

Struct to store Forward Sensitivity Finite State Projection outputs based on the sparse representation of the FSP solution.

# Fields 
- `t`: Array of solution output times.
- `p`: Array of `FspVectorSparse` instances. `p[i]` is the CME solution at time `t[i]`.
- 'S': Array of arrays of `FspVectorSparse` instances. `S[j][i]` is the `j`-th partial derivative of `p[i]`. 
- `sinks`: Probability mass accumulated at the sink states over time. `sinks[i]` is the vector of sink state probabilities at time `t[i]`.
- `dsinks`: Partial derivatives of sink probabilties. `dsinks[j][i]` is the `j`-th partial derivative of `sinks[i]`.

# Usage 
If `sol` is of type `FspOutputSparse`, `sol[i]` will return a slice, of type `FspOutputSliceSparse` of the solution set at the `i`-th index. 

# See also 
[`ForwardSensFspOutputSliceSparse`](@ref).

"""
Base.@kwdef mutable struct ForwardSensFspOutputSparse{NS,IntT<:Integer,RealT<:AbstractFloat}
    t::Vector{RealT}
    p::Vector{FspVectorSparse{NS,IntT,RealT}}
    sinks::Vector{AbstractVector{RealT}}
    S::Vector{Vector{FspVectorSparse{NS,IntT,RealT}}}
    dsinks::Vector{Vector{AbstractVector{RealT}}}
end

"""
    ForwardSensFspOutputSliceSparse

Struct to store the Forward Sensitivity FSP solution at a single time. Fields: `t`, `p`, `S`, `sinks`, `dsinks`.
"""
mutable struct ForwardSensFspOutputSliceSparse{NS,IntT<:Integer,RealT<:AbstractFloat}
    t::RealT
    p::FspVectorSparse{NS,IntT,RealT}
    sinks::AbstractVector{RealT}
    S::Vector{FspVectorSparse{NS,IntT,RealT}}
    dsinks::Vector{AbstractVector{RealT}}
end

function Base.getindex(sfspoutput::ForwardSensFspOutputSparse, ind::Integer)
    if ind > length(sfspoutput.t)
        throw(ArgumentError("Requested index exceeds array limit."))
    end
    return ForwardSensFspOutputSliceSparse(
        sfspoutput.t[ind],
        sfspoutput.p[ind],
        sfspoutput.sinks[ind],
        sfspoutput.S[ind],
        sfspoutput.dsinks[ind]
    )
    return 
end

function Base.getindex(sensfspoutput::ForwardSensFspOutputSparse, inds::Vector{T}) where {T <: Integer}    
    return [sensfspoutput[i] for i in inds]    
end

function Base.lastindex(sensfspoutput::ForwardSensFspOutputSparse)
    return length(sensfspoutput.t)
end

function Base.length(sfspoutput::ForwardSensFspOutputSparse)
    return length(sfspoutput.t)
end
function Base.size(sfspoutput::ForwardSensFspOutputSparse)
    return size(sfspoutput.t)
end

using Printf 
function Base.show(io::IO, output::ForwardSensFspOutputSparse)
    println(io, "Forward Sensitivity FSP solution set with $(length(output)) snapshots")
end

function Base.show(io::IO, outputslice::ForwardSensFspOutputSliceSparse)
    println(io, "Forward Sensitivity FSP solution output slice:")
    @printf(io, "t = %.2e \n", outputslice.t)
    println(io, "sinks = $(outputslice.sinks)")
    @printf(io, "p = \n")
    show(io, outputslice.p)    
    @printf(io, "S = \n")
    for svec in outputslice.S 
        show(io, svec)
    end
end