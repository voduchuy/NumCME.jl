export SparseFspOutput, SparseFspOutputSlice

Base.@kwdef struct SparseFspOutput{NS,IntT<:Integer,RealT<:AbstractFloat}
    t::Vector{RealT}
    p::Vector{SparseMultIdxVector{NS,IntT,RealT}}
    sinks::Vector{Vector{RealT}}
end

struct SparseFspOutputSlice
    t::AbstractFloat 
    p::SparseMultIdxVector
    sinks::Vector
end

function SparseFspOutput{NS}() where {NS, IntT <: Integer, RealT <: AbstractFloat}
    return SparseFspOutput{NS, Int64, Float64}(
        t = Vector{Float64}(),
        p = Vector{SparseMultIdxVector{NS, IntT, RealT}}(),
        sinks = Vector{Vector{RealT}}()
    )
end

function Base.getindex(fspoutput::SparseFspOutput, ind::Integer)
    if ind > length(fspoutput.t)
        throw(ArgumentError("Requested index exceeds array limit."))
    end
    return SparseFspOutputSlice(
        fspoutput.t[ind],
        fspoutput.p[ind],
        fspoutput.sinks[ind]
    )
    return 
end

function Base.getindex(fspoutput::SparseFspOutput, inds::Vector{T}) where {T <: Integer}    
    return [fspoutput[i] for i in inds]    
end

function Base.lastindex(fspoutput::SparseFspOutput)
    return lastindex(fspoutput.t)
end

function Base.length(fspoutput::SparseFspOutput)
    return length(fspoutput.t)
end
function Base.size(fspoutput::SparseFspOutput)
    return size(fspoutput.t)
end

using Printf 
function Base.show(io::IO, outputslice::SparseFspOutputSlice)
    println(io, "FSP solution output slice:")
    @printf(io, "t = %.2e \n", outputslice.t)
    println(io, "sinks = $(outputslice.sinks)")
    @printf(io, "p = \n")
    show(outputslice.p)    
end