
export RStepAdapter, SelectiveRStepAdapter

"""
Simple adapter based on reachability. Whenever the current FSP solution error is found to exceed the acceptable tolerance, the adapter will expand the state space by exploring all states that could be reached from the current state space within a set number of reaction events.

This adapter only works with `SparseStateSpace`.
"""
struct RStepAdapter <: SparseSpaceAdapter
    initial_step_count::Int
    max_step_count::Int
    dropstates::Bool
end

"""
`init!(statespace::SparseStateSpace, adapter::RStepAdapter, p::Vector{RealT}, t::RealT, fsptol::RealT) where {RealT <: AbstractFloat}`

Make the state space and solution vector ready for FSP integration.
"""
function init!(statespace::SparseStateSpace, adapter::RStepAdapter, p::Vector{RealT}, t::RealT, fsptol::RealT) where {RealT<:AbstractFloat}
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
    fsptol::RealT; integrator::Union{DEIntegrator,Nothing} = nothing) where {RealT<:AbstractFloat}

    if adapter.dropstates
        pids = sortperm(p)
        dropcount = sum((sum(p) .- cumsum(p[pids])) .>= (1.0 - t * fsptol / tend))
        dropids = sort(pids[1:dropcount])
        deleteat!(statespace, dropids)
        deleteat!(p, dropids)
    end

    nold = get_state_count(statespace)
    expand!(statespace, adapter.max_step_count)
    append!(p, zeros(RealT, get_state_count(statespace) - nold))
    nothing
end

"""
Adapter based on reachability. This adapter only explore new states through reaction channels with positive derivatives. 

This adapter only works with `SparseStateSpace`.
"""
struct SelectiveRStepAdapter <: SparseSpaceAdapter
    initial_step_count::Int
    max_step_count::Int
    dropstates::Bool   
end

"""
`init!(statespace::SparseStateSpace, adapter::SelectiveRStepAdapter, p::Vector{RealT}, t::RealT, fsptol::RealT) where {RealT <: AbstractFloat}`

Make the state space and solution vector ready for FSP integration.
"""
function init!(statespace::SparseStateSpace, adapter::SelectiveRStepAdapter, p::Vector{RealT}, t::RealT, fsptol::RealT) where {RealT<:AbstractFloat}
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
    integrator::Union{DEIntegrator,Nothing} = nothing) where {RealT<:AbstractFloat}
    if adapter.dropstates
        pids = sortperm(p)
        dropcount = sum((sum(p) .- cumsum(p[pids])) .> (1.0 - t * fsptol / tend))
        dropids = sort(pids[1:dropcount])
        deleteat!(statespace, dropids)
        deleteat!(p, dropids)        
    end
    sink_count = get_sink_count(statespace)
    du = similar(integrator.u)
    integrator.f(du, integrator.u, [], t)
    dsinks = du[end-sink_count+1:end]    
    expandreactions = findall(dsinks .> 0)    
    nold = get_state_count(statespace)
    expand!(statespace, adapter.max_step_count; onlyreactions = expandreactions)
    append!(p, zeros(RealT, get_state_count(statespace) - nold))               
    nothing
end