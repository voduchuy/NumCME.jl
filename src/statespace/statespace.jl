module statespace

using DataStructures: Deque
using StaticArrays

export FspStateSpace, FspStateSpaceBasic, rstep_expand!, get_state_count, get_sink_count

"""
Abstract type for FSP state space. This is the supertype of all concrete FSP state space implementations.
"""
abstract type FspStateSpace{NS, NR, T<:Int} end

"""
Basic FSP State space.
"""
Base.@kwdef mutable struct FspStateSpaceBasic{NS, NR, T<:Int} <: FspStateSpace{NS, NR, T}
    "Stoichiometry matrix S = [s₁ ... sₘ] of size N x M where N is the number of species, M the number of the reactions."
    stoich_mat::Matrix{T}

    "Number of sinks"
    sink_count::T

    "Array of CME states included in the subspace"
    states::Vector{MVector{NS, T}}

    "Dictionary of states, containing pairs `(x=>i)` for (i,x) in enumerate(states). The implementation must ensure that each state in `states` is a key in `state2idx` and conversely every key in `state2idx` exists in `states`."
    state2idx::Dict{MVector{NS, T},T}

    "List of state connectivity information. `state_connectivity[i][k] = j` if xᵢ = xⱼ + sₖ, that is, `states[i] = states[j] + stoich_mat[:, k]`. If there is no existing state that can reach xᵢ via reaction k, the implementation must ensure that `state_connectivity[i][k] = 0`."
    state_connectivity::Vector{MVector{NR,T}}

    "Matrix to store reaction events by which the included states transit to outside of the projected state space"
    sink_connectivity::Vector{MVector{NR,T}}    
end


"""
`FspStateSpaceBasic(stoich_mat, states)`

Construct a basic FSP state space with stoichiometry matrix `stoich_mat` and initial list of states `initstates`.
"""
function FspStateSpaceBasic(stoich_mat::Matrix{T}, initstates::Vector{Vector{T}}) where {T<:Int}
    species_count = size(stoich_mat, 1)
    reaction_count = size(stoich_mat, 2)
    sink_count = size(stoich_mat, 2)
    
    states = Vector{MVector{species_count, T}}()
    state2idx = Dict{MVector{species_count, T}, T}()    
    state_connectivity = Vector{MVector{reaction_count, T}}()    
    sink_connectivity = Vector{MVector{reaction_count, T}}()

    fspstatespace = FspStateSpaceBasic{species_count, reaction_count, T}(
        stoich_mat = stoich_mat,
        states = states,
        state2idx = state2idx,
        sink_count = sink_count,
        state_connectivity = state_connectivity,
        sink_connectivity = sink_connectivity,
    )

    for eachstate in initstates 
        _addstate!(MVector{species_count, T}(eachstate), fspstatespace)
    end    

    return fspstatespace
end

"""
`FspStateSpaceBasic(stoich_mat, init_state)` 

Construct a basic FSP state space with stoichiometry matrix `stoich_mat` and a single state `init_state`.
"""
function FspStateSpaceBasic(stoich_mat::Matrix{T}, init_state::Vector{T}) where {T<:Int}
    return FspStateSpaceBasic(stoich_mat, [init_state])
end

"""
Get number of states.
"""
function get_state_count(statespace::FspStateSpace)
    return length(statespace.states)
end

"""
Get number of sinks. 
"""
function get_sink_count(statespace::FspStateSpace)
    return statespace.sink_count
end

"""
`rstep_expand!(stateset::FspStateSpaceBasic, num_reachable_steps)`

Expand the FSP state space to include all states that are reachable from the existing states in `r` or fewer reaction events, where `r == num_reachable_steps`.
"""
function rstep_expand!(statespace::FspStateSpaceBasic{NS, NR, T}, num_reachable_steps::Int) where {NS, NR, T<:Int}
    if num_reachable_steps <= 0
        return 
    end    
    # Reference to mutable fields within stateset
    stoich_mat = statespace.stoich_mat
    states = statespace.states # Array of states 
    state2idx = statespace.state2idx # Dictionary of states
    sink_connectivity = statespace.sink_connectivity

    reaction_count = size(stoich_mat, 2)
    
    explorables = Deque{MVector{2,T}}() # State exploration queue    
    for idx = 1:length(states)
        for ir = 1:reaction_count
            if sink_connectivity[idx][ir] != 0
                push!(explorables, [idx, 0])
                break
            end
        end
    end 

    while !isempty(explorables)
        # Pop out an explorable state to explore         
        (idx, level) = popfirst!(explorables)
        
        (level < num_reachable_steps) && for stoichvec in eachcol(stoich_mat)
            candidate = states[idx] .+ stoichvec
            if is_all_nonnegative_(candidate)
                # Look up candidate state in the state dictionary. If the candidate does not exist in the dictionary, proceed to add it to the state list and update the parent state's connectivity.
                if get(state2idx, candidate, 0) == 0 
                    _addstate!(MVector{NS, T}(candidate), statespace)                
                    push!(explorables, [length(states), level+1])
                end
            end
        end
    end
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
function _addstate!(newstate::MVector{NS,T},
    statespace::FspStateSpaceBasic{NS, NR, T}) where {NS, NR, T<:Int}

    states = statespace.states
    state2idx = statespace.state2idx
    num_sinks = statespace.sink_count
    stoich_mat = statespace.stoich_mat
    state_connectivity = statespace.state_connectivity
    sink_connectivity = statespace.sink_connectivity    
    
    reaction_count = size(stoich_mat, 2)
    if get(state2idx, newstate, 0) == 0
        newidx = length(states)+1        
        push!(states, newstate)
        state2idx[newstate] = newidx

        # Determine connectivity between states         
        push!(state_connectivity, zeros(T, reaction_count))
        for ir = 1:reaction_count
            reachablestate = newstate - stoich_mat[:, ir]
            ridx = get(state2idx, reachablestate, 0)
            if ridx ≠ 0
                state_connectivity[newidx][ir] = ridx 
                sink_connectivity[ridx][ir] = 0
            end
        end

        # Determine states that can transit outside of the current truncated state space                
        isexplorable = false
        push!(sink_connectivity, zeros(Bool, num_sinks))
        for ir = 1:reaction_count
            reachablestate = newstate + stoich_mat[:, ir]
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

end