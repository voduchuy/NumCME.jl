export AbstractMultIdxVector, SparseMultIdxVector, sum

abstract type AbstractMultIdxVector end 

Base.@kwdef mutable struct SparseMultIdxVector{NS, IntT<:Integer, RealT<:AbstractFloat} <: AbstractMultIdxVector 
    states::Vector{MVector{NS, IntT}}
    values::Vector{RealT}
    state2idx::Dict{MVector{NS, IntT}, IntT}
end

function SparseMultIdxVector(states::Vector{MVector{NS, IntT}}, values::Vector{RealT};checksizes::Bool=true) where {NS, IntT<:Integer, RealT <: AbstractFloat}
    checksizes && (length(states) ≠ length(values)) && throw(ArgumentError("State and value lists must have equal lengths."))
    return SparseMultIdxVector{NS, IntT, RealT}(
        states=deepcopy(states),
        values=deepcopy(values),
        state2idx=Dict([x=>i for (i,x) in enumerate(states)])
    )
end

function SparseMultIdxVector(statespace::AbstractSparseStateSpace{NS, NR, IntT}, statevalpairs::Vector{Pair{VecT, RealT}}) where {NS, NR, IntT <: Integer, VecT <: AbstractVector, RealT <: AbstractFloat}
    states = deepcopy(statespace.states) 
    state2idx = deepcopy(statespace.state2idx)
    values = zeros(RealT, get_state_count(statespace))
    for (x,v) in statevalpairs
        idx = get(state2idx, x, 0)
        (idx ≠ 0) && (values[idx] = v )
    end
    return SparseMultIdxVector{NS, IntT, RealT}(
        states=states,
        values=values,
        state2idx=state2idx 
    )
end

function Base.sum(p::SparseMultIdxVector)
    Base.sum(p.values)
end

using Printf 
function Base.show(io::IO, vec::SparseMultIdxVector)
    nz = length(vec.values)
    if nz == 0
        println(io, "Sparse multi-indexed vector with no nonzero entries.")
        return nothing 
    end
    println("Sparse multi-indexed vector with $(nz) nonzero entries:")
    for i in 1:min(nz, 5)
        print(io, "$(vec.states[i]) ↦ ")
        @printf("%.2e \n", vec.values[i])        
    end    
    (nz > 10) && println(repeat(" ", length(vec.states[1])+5), "...")
    for i in max(nz-5, 6):nz 
        print(io, "$(vec.states[i]) ↦ ")
        @printf("%.2e \n", vec.values[i])        
    end
    nothing 
end

