export FspSolveOutput, FspSolveOutputSlice

Base.@kwdef struct FspSolveOutput{NS,IntT<:Integer,RealT<:AbstractFloat}
    t::Vector{RealT}
    p::Vector{FspSparseVector{NS,IntT,RealT}}
    sinks::Vector{Vector{RealT}}

    function FspSolveOutput{NS, IntT, RealT}(;t = Vector{RealT}(),
        p = Vector{FspSparseVector{NS, IntT, RealT}}(),
        sinks = Vector{Vector{RealT}}()) where {NS, IntT <: Integer, RealT <: AbstractFloat}    
            new(t, p, sinks)        
    end
end

struct FspSolveOutputSlice
    t::AbstractFloat 
    p::FspSparseVector
    sinks::Vector
end

function FspSolveOutput{NS}() where {NS, IntT <: Integer, RealT <: AbstractFloat}
    return FspSolveOutput{NS, Int64, Float64}(
        t = Vector{Float64}(),
        p = Vector{FspSparseVector{NS, IntT, RealT}}(),
        sinks = Vector{Vector{RealT}}()
    )
end

function Base.getindex(fspoutput::FspSolveOutput, ind::Integer)
    if ind > length(fspoutput.t)
        throw(ArgumentError("Requested index exceeds array limit."))
    end
    return FspSolveOutputSlice(
        fspoutput.t[ind],
        fspoutput.p[ind],
        fspoutput.sinks[ind]
    )
    return 
end

function Base.getindex(fspoutput::FspSolveOutput, inds::Vector{T}) where {T <: Integer}    
    return [fspoutput[i] for i in inds]    
end

function Base.length(fspoutput::FspSolveOutput)
    return length(fspoutput.t)
end
function Base.size(fspoutput::FspSolveOutput)
    return size(fspoutput.t)
end