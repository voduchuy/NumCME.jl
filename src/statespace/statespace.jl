using DataStructures: Deque

export AbstractStateSpace, AbstractSparseStateSpace, SparseStateSpace, expand!, get_state_count, get_sink_count

"""
Abstract type for FSP state space. This is the supertype of all concrete FSP state space implementations.
"""
abstract type AbstractStateSpace end
abstract type AbstractSparseStateSpace{NS,NR,IntT<:Integer} <: AbstractStateSpace end

"""
Basic Sparse FSP State space.
"""
Base.@kwdef mutable struct SparseStateSpace{NS,NR,IntT<:Integer} <: AbstractSparseStateSpace{NS,NR,IntT}
    "Stoichiometry matrix S = [s₁ ... sₘ] of size N x M where N is the number of species, M the number of the reactions."
    stoich_matrix::Matrix{IntT}

    "Number of sinks"
    sink_count::IntT

    "Array of CME states included in the subspace"
    states::Vector{MVector{NS,IntT}}

    "Dictionary of states, containing pairs `(x=>i)` for (i,x) in enumerate(states). The implementation must ensure that each state in `states` is a key in `state2idx` and conversely every key in `state2idx` exists in `states`."
    state2idx::Dict{MVector{NS,IntT},IntT}

    "List of state connectivity information. `state_connectivity[i][k] = j` if xᵢ = xⱼ + sₖ, that is, `states[i] = states[j] + stoich_mat[:, k]`. If there is no existing state that can reach xᵢ via reaction k, the implementation must ensure that `state_connectivity[i][k] = 0`."
    state_connectivity::Vector{MVector{NR,IntT}}

    "Matrix to store reaction events by which the included states transit to outside of the projected state space"
    sink_connectivity::Vector{MVector{NR,IntT}}
end
"""
`get_stoich_matrix(space::SparseStateSpace)`

Return the stoichiometry matrix.
"""
get_stoich_matrix(space::SparseStateSpace) = space.stoich_matrix

"""
`get_state_count(statespace::AbstractStateSpace)`

Return number of states.
"""
get_state_count(statespace::AbstractStateSpace) = length(statespace.states)

"""
`get_sink_count(statespace::AbstractStateSpace)`

Return number of sinks. 
"""
get_sink_count(statespace::AbstractStateSpace) = statespace.sink_count

"""
`get_states(statespace::AbstractSparseStateSpace)`

Return  
"""
get_states(statespace::AbstractStateSpace) = statespace.states 


"""
`SparseStateSpace(stoich_mat, states)`

Construct a basic FSP state space with stoichiometry matrix `stoich_mat` and initial list of states `initstates`.
"""
function SparseStateSpace(stoich_matrix::Matrix{IntT}, initstates::Vector) where {IntT<:Integer}
    species_count = size(stoich_matrix, 1)
    reaction_count = size(stoich_matrix, 2)
    sink_count = size(stoich_matrix, 2)

    states = Vector{MVector{species_count,IntT}}()
    state2idx = Dict{MVector{species_count,IntT},IntT}()
    state_connectivity = Vector{MVector{reaction_count,IntT}}()
    sink_connectivity = Vector{MVector{reaction_count,IntT}}()

    fspstatespace = SparseStateSpace{species_count,reaction_count,IntT}(
        stoich_matrix = stoich_matrix,
        states = states,
        state2idx = state2idx,
        sink_count = sink_count,
        state_connectivity = state_connectivity,
        sink_connectivity = sink_connectivity,
    )

    for eachstate in initstates
        _addstates!(fspstatespace, MVector{species_count,IntT}(eachstate))
    end

    return fspstatespace
end

"""
`SparseStateSpace(stoich_mat, init_state)` 

Construct a basic FSP state space with stoichiometry matrix `stoich_mat` and a single state `init_state`.
"""
function SparseStateSpace(stoich_matrix::Matrix{IntT}, init_state::Vector{IntT}) where {IntT<:Integer}
    return SparseStateSpace(stoich_matrix, [init_state])
end



"""
Expand the FSP state space to include all states that are reachable from the existing states in `r` or fewer reaction events, where `r == num_reachable_steps`.
"""
function expand!(statespace::SparseStateSpace{NS,NR,IntT}, expansionlevel::Integer; onlyreactions = []) where {NS,NR,IntT<:Integer}
    if expansionlevel <= 0
        return
    end
    # Reference to mutable fields within stateset
    stoich_matrix = statespace.stoich_matrix
    states = statespace.states # Array of states 
    state2idx = statespace.state2idx # Dictionary of states
    sink_connectivity = statespace.sink_connectivity

    reaction_count = size(stoich_matrix, 2)
    expandreactions = isempty(onlyreactions) ? (1:reaction_count) : onlyreactions

    explorables = Deque{MVector{2,IntT}}() # State exploration queue    
    for idx in 1:length(states)
        for ir in expandreactions
            if sink_connectivity[idx][ir] != 0
                push!(explorables, [idx, 0])
                break
            end
        end
    end

    while !isempty(explorables)
        # Pop out an explorable state to explore         
        (idx, level) = popfirst!(explorables)

        (level < expansionlevel) && for stoichvec in eachcol(stoich_matrix[:, expandreactions])
            candidate = states[idx] .+ stoichvec
            if is_all_nonnegative_(candidate)
                # Look up candidate state in the state dictionary. If the candidate does not exist in the dictionary, proceed to add it to the state list and update the parent state's connectivity.
                if get(state2idx, candidate, 0) == 0
                    _addstates!(statespace, MVector{NS,IntT}(candidate))
                    push!(explorables, [length(states), level + 1])
                end
            end
        end
    end
    nothing
end

"""
`is_all_positive_(x)`

Helper function. Returns `true` if all elements of a vector is positive. Otherwise returns `false`.
"""
function is_all_nonnegative_(x::AbstractVector)::Bool
    for xi in x
        xi < 0 && return false
    end
    return true
end

"""
Helper function to add a state vector `newstate` to the current state space `statespace`.
If `newstate` is already included in `statespace`, this function will return without modifying any of the input arguments.
"""
function _addstates!(statespace::SparseStateSpace{NS,NR,IntT}, newstate::MVector{NS,IntT}) where {NS,NR,IntT<:Integer}

    states = statespace.states
    state2idx = statespace.state2idx
    num_sinks = statespace.sink_count
    stoich_matrix = statespace.stoich_matrix
    state_connectivity = statespace.state_connectivity
    sink_connectivity = statespace.sink_connectivity

    reaction_count = size(stoich_matrix, 2)
    if get(state2idx, newstate, 0) == 0
        newidx = length(states) + 1
        push!(states, newstate)
        state2idx[newstate] = newidx

        # Determine connectivity between states         
        push!(state_connectivity, zeros(IntT, reaction_count))
        for ir = 1:reaction_count
            reachablestate = newstate - stoich_matrix[:, ir]
            ridx = get(state2idx, reachablestate, 0)
            if ridx ≠ 0
                state_connectivity[newidx][ir] = ridx
                sink_connectivity[ridx][ir] = 0
            end
        end

        # Determine states that can transit outside of the current truncated state space                
        push!(sink_connectivity, zeros(Bool, num_sinks))
        for ir = 1:reaction_count
            reachablestate = newstate + stoich_matrix[:, ir]
            if !(is_all_nonnegative_(reachablestate))
                continue
            end
            ridx = get(state2idx, reachablestate, 0)
            if ridx == 0
                sink_connectivity[newidx][ir] = ir
            else
                state_connectivity[ridx][ir] = newidx
            end
        end
    end
end
