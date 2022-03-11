using DataStructures: Deque 

export AbstractStateSpaceSparse, StateSpaceSparse, expand!, deleteat!, get_state_count, get_sink_count, get_stoich_matrix, get_statedict, get_state_connectivity, get_sink_connectivity

abstract type AbstractStateSpaceSparse{NS,NR,IntT<:Integer,SizeT<:Integer} <: AbstractStateSpace end

"""
    StateSpaceSparse{NS,NR,IntT<:Integer,SizeT<:Integer} <: AbstractStateSpaceSparse{NS,NR,IntT,SizeT}

Basic Sparse FSP State space.

# Fields 

    stoich_matrix::Matrix{IntT}
Stoichiometry matrix S = [s₁ ... sₘ] of size N x M where N is the number of species, M the number of the reactions.

    sink_count::SizeT
Number of sinks.

    states::Vector{MVector{NS,IntT}}
Array of CME states included in the subspace.

    state2idx::Dict{MVector{NS,IntT},SizeT}
Dictionary of states, containing pairs `(x=>i)` for (i,x) in enumerate(states). The implementation must ensure that each state in `states` is a key in `state2idx` and conversely every key in `state2idx` exists in `states`.

    state_connectivity::Vector{MVector{NR,SizeT}}
List of state connectivity information. `state_connectivity[i][k] = j` if xᵢ = xⱼ + sₖ, that is, `states[i] = states[j] + stoich_mat[:, k]`. If there is no existing state that can reach xᵢ via reaction k, the implementation must ensure that `state_connectivity[i][k] = 0`.

    sink_connectivity::Vector{MVector{NR,SizeT}}
Matrix to store reaction events by which the included states transit to outside of the projected state space.

# See also 
`expand!`,`deleteat!`
"""
Base.@kwdef mutable struct StateSpaceSparse{NS,NR,IntT<:Integer,SizeT<:Integer} <: AbstractStateSpaceSparse{NS,NR,IntT,SizeT}    
    "Stoichiometry matrix S = [s₁ ... sₘ] of size N x M where N is the number of species, M the number of the reactions."
    stoich_matrix::Matrix{IntT}

    "Number of sinks"
    sink_count::SizeT

    "Array of CME states included in the subspace"
    states::Vector{MVector{NS,IntT}}

    "Dictionary of states, containing pairs `(x=>i)` for (i,x) in enumerate(states). The implementation must ensure that each state in `states` is a key in `state2idx` and conversely every key in `state2idx` exists in `states`."
    state2idx::Dict{MVector{NS,IntT},SizeT}

    "List of state connectivity information. `state_connectivity[i][k] = j` if xᵢ = xⱼ + sₖ, that is, `states[i] = states[j] + stoich_mat[:, k]`. If there is no existing state that can reach xᵢ via reaction k, the implementation must ensure that `state_connectivity[i][k] = 0`."
    state_connectivity::Vector{MVector{NR,SizeT}}

    "Matrix to store reaction events by which the included states transit to outside of the projected state space"
    sink_connectivity::Vector{MVector{NR,SizeT}}
end

# Getters 
"""
`get_stoich_matrix(space::StateSpaceSparse)`

Return the stoichiometry matrix.
"""
get_stoich_matrix(space::StateSpaceSparse) = space.stoich_matrix

"""
`get_state_count(statespace::AbstractStateSpace)`

Return number of states.
"""
get_state_count(statespace::StateSpaceSparse) = length(statespace.states)

"""
`get_sink_count(statespace::AbstractStateSpace)`

Return number of sinks. 
"""
get_sink_count(statespace::StateSpaceSparse) = statespace.sink_count

"""
`get_states(statespace::AbstractStateSpaceSparse)`

Return list of states.
"""
get_states(statespace::StateSpaceSparse) = statespace.states

"""
`get_statedict(statespace::StateSpaceSparse)`

Return state dictionary.
"""
get_statedict(statespace::StateSpaceSparse) = statespace.state2idx

get_state_connectivity(statespace::StateSpaceSparse) = statespace.state_connectivity

get_sink_connectivity(statespace::StateSpaceSparse) = statespace.sink_connectivity

# Constructors 
"""
`StateSpaceSparse(stoich_mat, states; index_type::Type{<:Integer}=UInt32)`

Construct a basic FSP state space with stoichiometry matrix `stoich_mat` and initial list of states `initstates` with their integer entries being stored in type `IntT <: Integer`.

# Examples

```jldoctest
julia> S = [[1,0] [-1,0] [0,1] [0,-1]]
2×4 Matrix{Int64}:
 1  -1  0   0
 0   0  1  -1
 julia> states = [[0,1], [10, 1], [0, 10]]
 3-element Vector{Vector{Int64}}:
 [0, 1]
 [10, 1]
 [0, 10]
julia> StateSpaceSparse(S, states)
```

"""
function StateSpaceSparse(stoich_matrix::Matrix{IntT}, initstates::Vector; index_type::Type{SizeT} = UInt32) where {IntT<:Integer,SizeT<:Integer}
    species_count = size(stoich_matrix, 1)
    reaction_count = size(stoich_matrix, 2)
    sink_count = size(stoich_matrix, 2)

    states = Vector{MVector{species_count,IntT}}()
    state2idx = Dict{MVector{species_count,IntT},IntT}()
    state_connectivity = Vector{MVector{reaction_count,IntT}}()
    sink_connectivity = Vector{MVector{reaction_count,IntT}}()

    fspstatespace = StateSpaceSparse{species_count,reaction_count,IntT,index_type}(
        stoich_matrix = stoich_matrix,
        states = states,
        state2idx = state2idx,
        sink_count = sink_count,
        state_connectivity = state_connectivity,
        sink_connectivity = sink_connectivity,
    )

    _addstates!(fspstatespace, initstates)
    fspstatespace
end

"""
`StateSpaceSparse(stoich_mat::Matrix{IntT}, init_state::Vector{IntT}; index_type::Type{<:Integer}=UInt32)` 

Construct a basic FSP state space with stoichiometry matrix `stoich_mat` and a single state `init_state` with their integer entries being stored in type `IntT <: Integer`. The optional keyword argument `index_type` allows for more customization on internal indexing representations.

# Examples

```jldoctest
julia> S = [[1,0] [-1,0] [0,1] [0,-1]]
2×4 Matrix{Int64}:
 1  -1  0   0
 0   0  1  -1
 julia> x0 = [0,1]
julia> StateSpaceSparse(S, states)
"""
function StateSpaceSparse(stoich_matrix::Matrix{IntT}, init_state::Vector{IntT}; index_type::Type{SizeT} = UInt32) where {IntT<:Integer,SizeT<:Integer}
    return StateSpaceSparse(stoich_matrix, [init_state])
end



"""
Expand the FSP state space to include all states that are reachable from the existing states in `r` or fewer reaction events, where `r == num_reachable_steps`. The optional keyword argument `index_type` allows for more customization on internal indexing representations.
"""
function expand!(statespace::StateSpaceSparse{NS,NR,IntT,SizeT}, expansionlevel::Integer; onlyreactions = []) where {NS,NR,IntT<:Integer,SizeT<:Integer}
    if expansionlevel <= 0
        return
    end
    # Reference to mutable fields within stateset
    stoich_matrix = get_stoich_matrix(statespace)
    states = get_states(statespace)
    sink_connectivity = get_sink_connectivity(statespace)

    reaction_count = size(stoich_matrix, 2)
    expandreactions = isempty(onlyreactions) ? (1:reaction_count) : onlyreactions

    explorables = Deque{SizeT}()
    for idx in 1:length(states)
        for ir in expandreactions
            if sink_connectivity[idx][ir] != 0
                push!(explorables, idx)
                break
            end
        end
    end

    candidate_states = Vector{MVector{NS,IntT}}()
    for level in 1:expansionlevel
        resize!(candidate_states, length(explorables) * length(expandreactions))
        i = 1
        while !isempty(explorables)
            idx = pop!(explorables)
            @simd for ir in expandreactions
                candidate_states[i] = states[idx] .+ stoich_matrix[:, ir]
                i += 1
            end
        end
        statecount_old = length(states)
        _addstates!(statespace, candidate_states)
        statecount = length(states)
        for idx in statecount_old+1:statecount
            push!(explorables, idx)
        end
    end
    nothing
end

"""
`is_all_positive_(x)`

Helper function. Returns `true` if all elements of a vector is positive. Otherwise returns `false`.
"""
function _is_all_nonnegative(x::AbstractVector)::Bool
    for xi in x
        xi < 0 && return false
    end
    true
end

"""
Helper function to add a state vector `newstate` to the current state space `statespace`.
If `newstate` is already included in `statespace`, this function will return without modifying any of the input arguments.
"""
function _addstates!(statespace::StateSpaceSparse{NS,NR,IntT,SizeT}, newstates::Vector{VT}) where {NS,NR,IntT<:Integer,VT<:Union{Vector{IntT},MVector{NS,IntT}},SizeT<:Integer}
    states = get_states(statespace)
    state2idx = get_statedict(statespace)
    num_sinks = get_sink_count(statespace)
    stoich_matrix = get_stoich_matrix(statespace)
    state_connectivity = get_state_connectivity(statespace)
    sink_connectivity = get_sink_connectivity(statespace)
    reaction_count = size(stoich_matrix, 2)

    unique_newstates = Vector{MVector{NS,IntT}}()
    sizehint!(unique_newstates, length(newstates))
    newidx = length(states)
    for state in newstates
        if get(state2idx, state, 0) == 0 && _is_all_nonnegative(state)
            newidx += 1
            push!(unique_newstates, state)
            state2idx[state] = newidx
        end
    end

    statecount_old = length(states)
    added_state_count = length(unique_newstates)
    append!(states, unique_newstates)
    resize!(state_connectivity, statecount_old + added_state_count)
    resize!(sink_connectivity, statecount_old + added_state_count)
    statecount = length(states)

    @simd for newidx in statecount_old+1:statecount
        @inbounds state_connectivity[newidx] = zeros(IntT, reaction_count)
        @inbounds sink_connectivity[newidx] = zeros(IntT, reaction_count)
    end

    reachablestate = similar(states[1])
    @simd for newidx in statecount_old+1:statecount
        # Determine connectivity between states                 
        for ir = 1:reaction_count
            copy!(reachablestate, states[newidx])
            axpy!(-1, stoich_matrix[:, ir], reachablestate)
            # reachablestate = states[newidx] - stoich_matrix[:, ir]
            ridx = get(state2idx, reachablestate, 0)
            if ridx ≠ 0
                @inbounds state_connectivity[newidx][ir] = ridx
                @inbounds sink_connectivity[ridx][ir] = 0
            end
        end

        # Determine states that can transit outside of the current truncated state space                
        for ir = 1:reaction_count
            reachablestate = states[newidx] + stoich_matrix[:, ir]
            if _is_all_nonnegative(reachablestate)
                ridx = get(state2idx, reachablestate, 0)
                if ridx == 0
                    sink_connectivity[newidx][ir] = ir
                else
                    state_connectivity[ridx][ir] = newidx
                end
            end
        end
    end
end

import Base: deleteat!

"""
    deleteat!(statespace::StateSpaceSparse, ids::Vector{T}) where {T<:Integer}

Delete states with indices `ids` from the state space.
"""
function deleteat!(statespace::StateSpaceSparse, ids::Vector{T}) where {T<:Integer}
    stoich_matrix = get_stoich_matrix(statespace)
    reaction_count = size(stoich_matrix, 2)
    statecount_old = get_state_count(statespace)
    ids = unique(ids)

    # Map from current indices to new indices when the targeted states are removed
    newidxs = Vector{T}(1:statecount_old)
    newidxs[ids] .= 0
    idx = 1
    for i in 1:statecount_old
        if newidxs[i] ≠ 0
            newidxs[i] = idx
            idx += 1
        end
    end

    # Remove target states from member state list
    states = get_states(statespace)
    state2idx = get_statedict(statespace)
    state_connectivity = get_state_connectivity(statespace)
    sink_connectivity = get_sink_connectivity(statespace)

    for idx in ids
        delete!(state2idx, states[idx])
    end
    deleteat!(states, ids)
    deleteat!(state_connectivity, ids)
    deleteat!(sink_connectivity, ids)

    if length(states) == 0
        return nothing
    end
    # Reindexing
    reachablestate = similar(states[1])
    @simd for i in 1:length(states)
        @inbounds eachstate = states[i]
        state2idx[eachstate] = newidxs[state2idx[eachstate]]
        for ir in 1:reaction_count
            @inbounds state_connectivity[i][ir] = (oldidx = state_connectivity[i][ir]) ≠ 0 ? newidxs[oldidx] : 0
        end
    end
    @simd for i in 1:length(states)
        eachstate = states[i]
        for ir in 1:reaction_count
            if sink_connectivity[i][ir] == 0
                copyto!(reachablestate, eachstate)
                axpy!(1, stoich_matrix[:, ir], reachablestate)
                if _is_all_nonnegative(reachablestate) && (get(state2idx, reachablestate, 0) == 0)
                    @inbounds sink_connectivity[i][ir] = ir
                end
            end
        end
    end
    nothing
end
