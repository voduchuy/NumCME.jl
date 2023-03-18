import SparseArrays: nnz

export AbstractFspVector, FspVectorSparse, sum, get_states, get_values, get_statedict , nnz

abstract type AbstractFspVector end


"""
    FspVectorSparse{NS, IntT<:Integer, RealT<:AbstractFloat} <: AbstractFspVector 

Sparse multi-indexed vector (i.e., sparse "tensor") with entries of type `RealT` and multi-indices of type `IntT` for storing CME solution and sensitivity vectors output by the Finite State Projection.   
"""
Base.@kwdef mutable struct FspVectorSparse{NS,IntT<:Integer,RealT<:AbstractFloat} <: AbstractFspVector
    states::Vector{MVector{NS,IntT}}
    values::Vector{RealT}
    state2idx::Dict{MVector{NS,IntT},IntT}
end
get_states(v::FspVectorSparse) = v.states
get_values(v::FspVectorSparse) = v.values
get_statedict(v::FspVectorSparse) = v.state2idx
nnz(p::FspVectorSparse) = length(get_states(p))

"""
    FspVectorSparse(states::Vector{MVector{NS,IntT}}, values::Vector{RealT}; checksizes::Bool = true) where {NS,IntT<:Integer,RealT<:AbstractFloat}

Construct a sparse representation of a FSP-truncated probability distribution with support `states` and proability values `values`.
"""
function FspVectorSparse(states::Vector{MVector{NS,IntT}}, values::Vector{RealT}; checksizes::Bool = true) where {NS,IntT<:Integer,RealT<:AbstractFloat}
    checksizes && (length(states) ≠ length(values)) && throw(ArgumentError("State and value lists must have equal lengths."))
    return FspVectorSparse{NS,IntT,RealT}(
        states = deepcopy(states),
        values = deepcopy(values),
        state2idx = Dict([x => i for (i, x) in enumerate(states)])
    )
end

"""
    FspVectorSparse(statespace::AbstractStateSpaceSparse{NS,NR,IntT}, statevalpairs::Vector{Pair{VecT,RealT}}) where {NS,NR,IntT<:Integer,VecT<:AbstractVector,RealT<:AbstractFloat}

Construct a sparse representation of a FSP-truncated probability distribution on the state space `statespace` and proability values `values`.
"""
function FspVectorSparse(statespace::AbstractStateSpaceSparse{NS,NR,IntT}, statevalpairs::Vector{Pair{VecT,RealT}}) where {NS,NR,IntT<:Integer,VecT<:AbstractVector,RealT<:AbstractFloat}
    states = deepcopy(get_states(statespace))
    state2idx = deepcopy(get_statedict(statespace))
    values = zeros(RealT, get_state_count(statespace))
    for (x, v) in statevalpairs
        idx = get(state2idx, x, 0)
        (idx ≠ 0) && (values[idx] = v)
    end
    return FspVectorSparse{NS,IntT,RealT}(
        states = states,
        values = values,
        state2idx = state2idx
    )
end

function Base.sum(p::FspVectorSparse)
    Base.sum(p.values)
end

"""
    sum(p::FspVectorSparse, dims)

Reduce the FSP vector `p` by summing over species specified in `dims` and return a new FSP vector over a reduced-dimension state space.
"""
function Base.sum(p::FspVectorSparse{NS,IntT,RealT}, dims::AbstractVector{<:Integer})::FspVectorSparse where {NS,IntT<:Integer,RealT<:AbstractFloat}
    dims = unique(dims)
    !((min(dims...) >= 1) && (max(dims...) <= NS)) && throw(ArgumentError("Input dimensions must be between 1 and $(NS)."))

    keepdims = setdiff(1:NS, dims)

    full_states = get_states(p)
    full_vals = get_values(p)

    reduced_species_count = NS - length(dims)
    reduced_states = Vector{MVector{reduced_species_count,IntT}}()
    reduced_state2idx = Dict{MVector{reduced_species_count,IntT}, IntT}()
    reduced_vals = Vector{RealT}()

    rvec_length = 0    
    for (i, eachstate) in enumerate(full_states)
        reduced_state = eachstate[keepdims]
        reduced_idx = get(reduced_state2idx, reduced_state, 0)
        if reduced_idx == 0
            rvec_length += 1
            reduced_state2idx[reduced_state] = rvec_length
            push!(reduced_states, copy(reduced_state))
            push!(reduced_vals, full_vals[i])
        else
            reduced_vals[reduced_idx] += full_vals[i]
        end
    end

    return FspVectorSparse{reduced_species_count,IntT,RealT}(
        states = reduced_states,
        values = reduced_vals,
        state2idx = reduced_state2idx
    )
end


"""
    Array(p::FspVectorSparse)

Construct a new dense `N`-dimensional array from the FSP vector.
"""
function Base.Array(p::FspVectorSparse{NS,IntT,RealT}) where {NS, IntT<:Integer, RealT<:AbstractFloat}
    states = get_states(p)
    if isempty(states)
        throw(ArgumentError("Cannot construct dense array from empty FspVector instance."))
    end
    values = get_values(p)

    dimcount = length(states[1])
    bounds = zeros(Int, dimcount)
    for eachstate in states 
        for i in 1:dimcount 
            bounds[i] = max(bounds[i], eachstate[i]+1)
        end
    end

    out = zeros(RealT, bounds...)
    tmp = MVector{NS,IntT}([0 for i in 1:NS])
    for i in 1:length(values)
        tmp .= states[i] .+ 1
        out[tmp...] = values[i]
    end
    return out 
end

using Printf
function Base.show(io::IO, vec::FspVectorSparse)
    nz = length(vec.values)
    if nz == 0
        println(io, "Sparse multi-indexed vector with no nonzero entries.")
        return nothing
    end
    println("Sparse multi-indexed vector with $(nz) nonzero entries:")
    for i in 1:min(nz, 5)
        print(io, "$(vec.states[i]) ↦ ")
        @printf(io, "%.2e \n", vec.values[i])
    end
    (nz > 10) && println(repeat(" ", length(vec.states[1]) + 5), "...")
    for i in max(nz - 5, 6):nz
        print(io, "$(vec.states[i]) ↦ ")
        @printf(io, "%.2e \n", vec.values[i])
    end
    nothing
end


