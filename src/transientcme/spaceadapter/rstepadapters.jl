
export RStepAdapter, SelectiveRStepAdapter

"""
Simple adapter based on reachability. Whenever the current FSP solution error is found to exceed the acceptable tolerance, the adapter will expand the state space by exploring all states that could be reached from the current state space within a set number of reaction events.

This adapter only works with `SparseStateSpace`.
"""
struct RStepAdapter <: SparseSpaceAdapter
    initial_step_count::Int 
    max_step_count::Int 
end

"""
`init!(statespace::SparseStateSpace, adapter::RStepAdapter, p::Vector{RealT}, t::RealT, fsptol::RealT) where {RealT <: AbstractFloat}`

Make the state space and solution vector ready for FSP integration.
"""
function init!(statespace::SparseStateSpace, adapter::RStepAdapter, p::Vector{RealT}, t::RealT, fsptol::RealT) where {RealT <: AbstractFloat}
    nold = get_state_count(statespace)
    expand!(statespace, adapter.initial_step_count)
    append!(p, zeros(RealT, get_state_count(statespace) - nold))    
    nothing
end

"""
`adapt!(statespace::SparseStateSpace, adapter::RStepAdapter, p::Vector{RealT}, sinks::Vector{RealT}, t::RealT, fsptol::RealT; integrator::Union{DEIntegrator, Nothing}=nothing) where {RealT <: AbstractFloat}`

Adapt the state space and probability vector based on current error recorded in the `sinks` vector.
"""
function adapt!(statespace::SparseStateSpace, adapter::RStepAdapter, p::Vector{RealT}, 
    sinks::Vector{RealT}, 
    t::RealT, tend::RealT,
    fsptol::RealT; integrator::Union{DEIntegrator, Nothing}=nothing) where {RealT <: AbstractFloat}
    nold = get_state_count(statespace)
    expand!(statespace, adapter.max_step_count)
    append!(p, zeros(RealT, get_state_count(statespace) - nold))    
    nothing
end

"""
Adapter based on reachability. This adapter evaluates the probability mass that has flown 
through different reaction channels and only explore new states through reaction channels that accumulate the largest error. 

This adapter only works with `SparseStateSpace`.
"""
struct SelectiveRStepAdapter <: SparseSpaceAdapter
    initial_step_count::Int 
    max_step_count::Int
end

"""
`init!(statespace::SparseStateSpace, adapter::SelectiveRStepAdapter, p::Vector{RealT}, t::RealT, fsptol::RealT) where {RealT <: AbstractFloat}`

Make the state space and solution vector ready for FSP integration.
"""
function init!(statespace::SparseStateSpace, adapter::SelectiveRStepAdapter, p::Vector{RealT}, t::RealT, fsptol::RealT) where {RealT <: AbstractFloat}
    nold = get_state_count(statespace)
    expand!(statespace, adapter.initial_step_count)
    append!(p, zeros(RealT, get_state_count(statespace) - nold))    
    nothing
end

"""
`adapt!(statespace::SparseStateSpace, adapter::SelectiveRStepAdapter, p::Vector{RealT}, sinks::Vector{RealT}, t::RealT, fsptol::RealT; integrator::Union{DEIntegrator, Nothing}=nothing) where {RealT <: AbstractFloat}`

Adapt the state space and probability vector based on current error recorded in the `sinks` vector.
"""
function adapt!(statespace::SparseStateSpace, adapter::SelectiveRStepAdapter, 
    p::Vector{RealT}, sinks::Vector{RealT}, 
    t::RealT, tend::RealT, fsptol::RealT; 
    integrator::Union{DEIntegrator, Nothing}=nothing) where {RealT <: AbstractFloat}        
    expandreactions = findall(sinks .>= (t*fsptol/tend)/length(sinks))
    nold = get_state_count(statespace)
    expand!(statespace, adapter.max_step_count; onlyreactions=expandreactions)
    append!(p, zeros(RealT, get_state_count(statespace) - nold))    
    nothing
end