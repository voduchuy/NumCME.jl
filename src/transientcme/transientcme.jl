import DifferentialEquations
import Sundials
using DifferentialEquations.DiffEqBase: ODEProblem, AbstractODEAlgorithm
using SparseArrays

export TransientCmeAlgorithm, FixedSparseFsp, AdaptiveSparseFsp, solve

abstract type TransientCmeAlgorithm end

Base.@kwdef mutable struct FixedSparseFsp <: TransientCmeAlgorithm
    ode_method::Union{Nothing, AbstractODEAlgorithm}
end

Base.@kwdef mutable struct AdaptiveSparseFsp <: TransientCmeAlgorithm
    ode_method::Union{Nothing, AbstractODEAlgorithm}
    space_adapter::SparseSpaceAdapter
end

struct FspSolveOutput{NS,IntT<:Integer,RealT<:AbstractFloat}
    t::RealT
    p::FspSparseVector{NS,IntT,RealT}
    sinks::Vector{RealT}
end

function solve(model::CmeModel,
    initial_distribution::FspSparseVector,
    tspan::Union{Vector,Tuple},
    fspalgorithm::FixedSparseFsp;
    saveat=[],
    fsptol::AbstractFloat = 1.0E-6,
    odeatol::AbstractFloat = 1.0E-10,
    odertol::AbstractFloat = 1.0E-4)
    𝔛 = SparseStateSpace(model.stoich_matrix, initial_distribution.states)
    sink_count = get_sink_count(𝔛)
    p0 = initial_distribution.values

    A = FspSparseMatrix(𝔛, model.propensities)
    u0 = [p0; zeros(sink_count)]
    function odefun!(du, u, θ, t)
        matvec!(t, A, u, du)
        nothing
    end

    fspprob = ODEProblem(odefun!, u0, tspan)    
    solutions = Sundials.solve(fspprob, fspalgorithm.ode_method, atol = odeatol, rtol = odertol, saveat=saveat)

    outputs = []
    for (t, u) in zip(solutions.t, solutions.u)
        p = u[1:end-sink_count]
        sinks = u[end-sink_count+1:end]
        push!(
            outputs,
            FspSolveOutput(
                t,
                FspSparseVector(initial_distribution.states, p),
                sinks)
        )
    end

    return outputs
end

function solve(model::CmeModel, 
    initial_distribution::FspSparseVector, 
    tspan::Vector{AbstractFloat}, 
    fspalgorithm::AdaptiveSparseFsp; 
    saveat = [],
    fsptol::AbstractFloat = 1.0E-6, 
    odeatol::AbstractFloat = 1.0E-10, 
    odertol::AbstractFloat = 1.0E-4)
    
    tₙ = min(tspan)    
    tmax = max(tspan)
    𝔛 = SparseStateSpace(model.stoich_matrix, [e[1] for e in p0])
    sink_count = get_sink_count(𝔛)
        
    uₙ = [initial_distribution.values; zeros(sink_count)]
    outputs = []
    while tₙ < tmax
        A = FspSparseMatrix(𝔛, model.propensities)
        function odefun!(du, u, t, p)
            matvec!(t, A, u, du)
            nothing
        end

        fspprob = ODEProblem(odefun!, uₙ, tspan)    
        solutions = Sundials.solve(fspprob, fspalgorithm.ode_method, atol = odeatol, rtol = odertol, saveat=saveat)
            
    end
    return outputs 
end
end