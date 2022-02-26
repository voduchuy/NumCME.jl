import DifferentialEquations as DE
import Sundials
using DifferentialEquations.DiffEqBase: ODEProblem, AbstractODEAlgorithm
using SparseArrays

export TransientCmeAlgorithm, FixedSparseFsp, AdaptiveSparseFsp, solve

include("fspoutput.jl")
include("spaceadapter.jl")

abstract type TransientCmeAlgorithm end

Base.@kwdef mutable struct FixedSparseFsp <: TransientCmeAlgorithm
    ode_method::Union{Nothing,AbstractODEAlgorithm}
end

function solve(model::CmeModel,
    initial_distribution::FspSparseVector,
    tspan::Union{Vector,Tuple},
    fspalgorithm::FixedSparseFsp;
    saveat = [],
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
    solutions = DE.solve(fspprob, fspalgorithm.ode_method, atol = odeatol, rtol = odertol, saveat = saveat)

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

Base.@kwdef mutable struct AdaptiveSparseFsp <: TransientCmeAlgorithm
    ode_method::Union{Nothing,AbstractODEAlgorithm}
    space_adapter::SparseSpaceAdapter
end

function solve(model::CmeModel,
    initial_distribution::FspSparseVector{NS,IntT,RealT},
    tspan::Tuple{AbstractFloat,AbstractFloat},
    fspalgorithm::AdaptiveSparseFsp;
    saveat = [],
    fsptol::AbstractFloat = 1.0E-6,
    odeatol::AbstractFloat = 1.0E-10,
    odertol::AbstractFloat = 1.0E-4) where {NS,IntT<:Integer,RealT<:AbstractFloat}

    tstart = min(tspan...)
    tend = max(tspan...)
    adapter = fspalgorithm.space_adapter

    p0 = deepcopy(initial_distribution.values)
    𝔛 = SparseStateSpace(model.stoich_matrix, initial_distribution.states)
    sink_count = get_sink_count(𝔛)
    init!(𝔛, adapter, p0, tstart, fsptol)    

    tnow = tstart
    unow = [p0; zeros(sink_count)]
    A = FspSparseMatrix(𝔛, model.propensities)
    function fsprhs!(du, u, θ, t)
        matvec!(t, A, u, du)
        nothing
    end
    function fsp_error_constraint(u, t, integrator)
        sinks = u[end-sink_count:end]
        return sum(sinks) - fsptol * t / tend
    end
    function affect!(integrator)
        DE.terminate!(integrator)
    end
    fsperrorcb = DE.ContinuousCallback(
        fsp_error_constraint, affect!, save_positions = (false, false)
    )

    output = FspSolveOutput{NS,IntT,RealT}()
    while tnow < tend
        fspprob = DE.ODEProblem(fsprhs!, unow, (tnow, tend))

        integrator = DE.init(fspprob, fspalgorithm.ode_method, atol = odeatol, rtol = odertol, callback = fsperrorcb, saveat = saveat)

        DE.step!(integrator, tend - tnow, true)

        for (t, u) in zip(integrator.sol.t, integrator.sol.u)
            push!(output.t, t)
            push!(output.p, FspSparseVector(𝔛.states, u[1:end-sink_count]))
            push!(output.sinks, u[end-sink_count+1:end])
        end

        tnow = integrator.t
        if tnow < tend
            p = integrator.u[1:end-sink_count]
            sinks = integrator.u[end-sink_count+1:end]
            adapt!(𝔛, adapter, p, sinks, tnow, fsptol)
            unow = [p; sinks]
            A = FspSparseMatrix(𝔛, model.propensities)
        else
            u = integrator.u
            push!(output.t, tnow)
            push!(output.p, FspSparseVector(𝔛.states, u[1:end-sink_count]))
            push!(output.sinks, u[end-sink_count+1:end])
        end
    end

    return output
end
