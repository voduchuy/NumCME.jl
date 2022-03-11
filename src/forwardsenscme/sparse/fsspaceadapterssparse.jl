export ForwardSensRStepAdapter, AbstractForwardSensSpaceAdapterSparse

abstract type AbstractForwardSensSpaceAdapterSparse <: AbstractForwardSensSpaceAdapter end 

"""
Simple adapter based on reachability. Whenever the current FSP solution error is found to exceed the acceptable tolerance, the adapter will expand the state space by exploring all states that could be reached from the current state space within a set number of reaction events.

This adapter only works with `StateSpaceSparse`.
"""
struct ForwardSensRStepAdapter <: AbstractForwardSensSpaceAdapterSparse
    initial_step_count::Int
    max_step_count::Int
    dropstates::Bool
end

"""
`init!(statespace::StateSpaceSparse, adapter::RStepAdapter, p::Vector{RealT}, S::Vector{Vector{RealT}}, t::RealT, fsptol::RealT) where {RealT <: AbstractFloat}`

Initialize the state space, solution vector and sensitivity indices for FSP integration.
"""
function init!(statespace::StateSpaceSparse, adapter::ForwardSensRStepAdapter, p::Vector{RealT}, S::Vector{Vector{RealT}}, t::AbstractFloat, fsptol::AbstractFloat) where {RealT <: AbstractFloat}
    nold = get_state_count(statespace)
    expand!(statespace, adapter.initial_step_count)
    append!(p, zeros(RealT, get_state_count(statespace) - nold))
    for svec in S 
        append!(svec, zeros(RealT, get_state_count(statespace) - nold))
    end
    nothing
end

"""
`adapt!(statespace::StateSpaceSparse, adapter::RStepAdapter, p::Vector{RealT}, sinks::Vector{RealT}, t::RealT, fsptol::RealT; integrator::Union{DEIntegrator, Nothing}=nothing) where {RealT <: AbstractFloat}`

Adapt the state space, probability vector and sensitivity vectors based on current error recorded in the `sinks` vector.
"""
function adapt!(statespace::StateSpaceSparse, adapter::ForwardSensRStepAdapter, p::Vector{RealT}, S::Vector{Vector{RealT}},
    sinks::Vector{RealT},
    dsinks::Vector{Vector{RealT}},
    t::RealT, tend::RealT,
    fsptol::RealT; integrator::Union{DEIntegrator,Nothing} = nothing) where {RealT<:AbstractFloat}

    if adapter.dropstates
        pids = sortperm(p)
        dropcount = sum((sum(p) .- cumsum(p[pids])) .>= (1.0 - t * fsptol / tend))
        dropids = sort(pids[1:dropcount])
        deleteat!(statespace, dropids)
        deleteat!(p, dropids)
        for svec in S 
            deleteat!(svec, dropids)
        end
    end

    nold = get_state_count(statespace)
    expand!(statespace, adapter.max_step_count)
    append!(p, zeros(RealT, get_state_count(statespace) - nold))
    for svec in S 
        append!(svec, zeros(RealT, get_state_count(statespace) - nold))
    end
    nothing
end



