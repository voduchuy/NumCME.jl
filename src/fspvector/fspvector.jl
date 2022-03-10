export AbstractMultIdxVector, MultIdxVectorSparse, sum, get_states, get_values, get_statedict 

abstract type AbstractMultIdxVector end

"""
    MultIdxVectorSparse{NS, IntT<:Integer, RealT<:AbstractFloat} <: AbstractMultIdxVector 

Sparse multi-indexed vector (i.e., sparse "tensor") with entries of type `RealT` and multi-indices of type `IntT`.   
"""
Base.@kwdef mutable struct MultIdxVectorSparse{NS,IntT<:Integer,RealT<:AbstractFloat} <: AbstractMultIdxVector
    states::Vector{MVector{NS,IntT}}
    values::Vector{RealT}
    state2idx::Dict{MVector{NS,IntT},IntT}
end
get_states(v::MultIdxVectorSparse) = v.states
get_values(v::MultIdxVectorSparse) = v.values
get_statedict(v::MultIdxVectorSparse) = v.state2idx

"""
    MultIdxVectorSparse(states::Vector{MVector{NS, IntT}}, values::Vector{RealT})

Construct a sparse representation of a FSP-truncated probability distribution with support `states` and proability values `values`.
"""
function MultIdxVectorSparse(states::Vector{MVector{NS,IntT}}, values::Vector{RealT}; checksizes::Bool = true) where {NS,IntT<:Integer,RealT<:AbstractFloat}
    checksizes && (length(states) ≠ length(values)) && throw(ArgumentError("State and value lists must have equal lengths."))
    return MultIdxVectorSparse{NS,IntT,RealT}(
        states = deepcopy(states),
        values = deepcopy(values),
        state2idx = Dict([x => i for (i, x) in enumerate(states)])
    )
end

"""
    MultIdxVectorSparse(statespace::AbstractStateSpaceSparse{NS, NR, IntT}, statevalpairs::Vector{Pair{VecT, RealT}})

Construct a sparse representation of a FSP-truncated probability distribution on the state space `statespace` and proability values `values`.
"""
function MultIdxVectorSparse(statespace::AbstractStateSpaceSparse{NS,NR,IntT}, statevalpairs::Vector{Pair{VecT,RealT}}) where {NS,NR,IntT<:Integer,VecT<:AbstractVector,RealT<:AbstractFloat}
    states = deepcopy(statespace.states)
    state2idx = deepcopy(statespace.state2idx)
    values = zeros(RealT, get_state_count(statespace))
    for (x, v) in statevalpairs
        idx = get(state2idx, x, 0)
        (idx ≠ 0) && (values[idx] = v)
    end
    return MultIdxVectorSparse{NS,IntT,RealT}(
        states = states,
        values = values,
        state2idx = state2idx
    )
end

function Base.sum(p::MultIdxVectorSparse)
    Base.sum(p.values)
end

function Base.sum(p::MultIdxVectorSparse{NS,IntT,RealT}, dims::AbstractVector{<:Integer})::MultIdxVectorSparse where {NS,IntT<:Integer,RealT<:AbstractFloat}
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
    reduced_state = MVector{reduced_species_count,IntT}(zeros(IntT, reduced_species_count))
    for (i, eachstate) in enumerate(full_states)
        reduced_state .= eachstate[keepdims]
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

    return MultIdxVectorSparse{reduced_species_count,IntT,RealT}(
        states = reduced_states,
        values = reduced_vals,
        state2idx = reduced_state2idx
    )
end

using Printf
function Base.show(io::IO, vec::MultIdxVectorSparse)
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


