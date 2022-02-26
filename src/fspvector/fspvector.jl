using SparseArrays
using StaticArrays

export AbstractFspVector, FspSparseVector, sum

abstract type AbstractFspVector end 

Base.@kwdef mutable struct FspSparseVector{NS, IntT<:Integer, RealT<:AbstractFloat} <: AbstractFspVector 
    states::Vector{MVector{NS, IntT}}
    values::Vector{RealT}
    state2idx::Dict{MVector{NS, IntT}, IntT}
end

function FspSparseVector(states::Vector{MVector{NS, IntT}}, values::Vector{RealT};checksizes::Bool=true) where {NS, IntT<:Integer, RealT <: AbstractFloat}
    checksizes && (length(states) ≠ length(values)) && throw(ArgumentError("State and value lists must have equal lengths."))
    return FspSparseVector{NS, IntT, RealT}(
        states=deepcopy(states),
        values=deepcopy(values),
        state2idx=Dict([x=>i for (i,x) in enumerate(states)])
    )
end

function FspSparseVector{RealT}(statespace::AbstractSparseStateSpace{NS, NR, IntT}, statevalpairs::Vector{Pair{VecT, RealT}}) where {NS, NR, IntT <: Integer, VecT <: AbstractVector, RealT <: AbstractFloat}
    states = deepcopy(statespace.states) 
    state2idx = deepcopy(statespace.state2idx)
    values = zeros(RealT, get_state_count(statespace))
    for (x,v) in statevalpairs
        idx = get(state2idx, x, 0)
        (idx ≠ 0) && (values[idx] = v )
    end
    return FspSparseVector{NS, IntT, RealT}(
        states=states,
        values=values,
        state2idx=state2idx 
    )
end

function FspSparseVector(statespace::AbstractSparseStateSpace{NS, NR, IntT}, statevalpairs::Vector{Pair{VecT, Float64}}) where {NS, NR, IntT <: Integer, VecT <: AbstractVector}
    return FspSparseVector{Float64}(statespace, statevalpairs)
end

function Base.sum(p::FspSparseVector)
    Base.sum(p.values)
end



