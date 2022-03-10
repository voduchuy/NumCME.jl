export FspOutputSparse, FspOutputSliceSparse

Base.@kwdef struct FspOutputSparse{NS,IntT<:Integer,RealT<:AbstractFloat}
    t::Vector{RealT}
    p::Vector{MultIdxVectorSparse{NS,IntT,RealT}}
    sinks::Vector{Vector{RealT}}
end

struct FspOutputSliceSparse
    t::AbstractFloat 
    p::MultIdxVectorSparse
    sinks::Vector
end

function FspOutputSparse{NS}() where {NS, IntT <: Integer, RealT <: AbstractFloat}
    return FspOutputSparse{NS, Int64, Float64}(
        t = Vector{Float64}(),
        p = Vector{MultIdxVectorSparse{NS, IntT, RealT}}(),
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