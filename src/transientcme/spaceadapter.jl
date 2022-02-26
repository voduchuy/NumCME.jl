export SparseSpaceAdapter, SpaceAdapter, RStepAdapter, SelectiveRStepAdapter, init!, adapt!

abstract type SpaceAdapter end
abstract type SparseSpaceAdapter <: SpaceAdapter end 
using DifferentialEquations.DiffEqBase: DEIntegrator

struct RStepAdapter <: SparseSpaceAdapter
    initial_step_count::Int 
    max_step_count::Int 
end

struct SelectiveRStepAdapter <: SparseSpaceAdapter
    initial_step_count::Int 
    max_step_count::Int
end

function init!(statespace::SparseStateSpace, adapter::RStepAdapter, p::Vector{RealT}, t::RealT, fsptol::RealT) where {RealT <: AbstractFloat}
    nold = get_state_count(statespace)
    expand!(statespace, adapter.initial_step_count)
    append!(p, zeros(RealT, get_state_count(statespace) - nold))    
    nothing
end

function adapt!(statespace::SparseStateSpace, adapter::RStepAdapter, p::Vector{RealT}, sinks::Vector{RealT}, t::RealT, fsptol::RealT; integrator::Union{DEIntegrator, Nothing}=nothing) where {RealT <: AbstractFloat}
    nold = get_state_count(statespace)
    expand!(statespace, adapter.max_step_count)
    append!(p, zeros(RealT, get_state_count(statespace) - nold))    
    nothing
end

function init!(statespace::SparseStateSpace, adapter::SelectiveRStepAdapter, p::Vector{RealT}, t::RealT, fsptol::RealT) where {RealT <: AbstractFloat}
    nold = get_state_count(statespace)
    expand!(statespace, adapter.initial_step_count)
    append!(p, zeros(RealT, get_state_count(statespace) - nold))    
    nothing
end

function adapt!(statespace::SparseStateSpace, adapter::SelectiveRStepAdapter, p::Vector{RealT}, sinks::Vector{RealT}, t::RealT, fsptol::RealT; integrator::Union{DEIntegrator, Nothing}=nothing) where {RealT <: AbstractFloat}
    expandreactions = findall(sinks .> sum(sinks)/length(sinks))
    nold = get_state_count(statespace)
    expand!(statespace, adapter.max_step_count; onlyreactions=expandreactions)
    append!(p, zeros(RealT, get_state_count(statespace) - nold))    
    nothing
end
