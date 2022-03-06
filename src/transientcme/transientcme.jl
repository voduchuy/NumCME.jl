import DifferentialEquations as DE
import Sundials
using DifferentialEquations.DiffEqBase: ODEProblem, AbstractODEAlgorithm


export TransientCmeAlgorithm, FixedSparseFsp, AdaptiveSparseFsp, solve

include("fspoutput.jl")
include("spaceadapter/spaceadapter.jl")

abstract type TransientCmeAlgorithm end

Base.@kwdef mutable struct FixedSparseFsp <: TransientCmeAlgorithm
    ode_method::Union{Nothing,AbstractODEAlgorithm}
end

function solve(model::CmeModel,
    initial_distribution::SparseMultIdxVector{NS,IntT,RealT},
    tspan::Union{Vector,Tuple},
    fspalgorithm::FixedSparseFsp; 
    saveat = [], fsptol::AbstractFloat = 1.0E-6, odeatol::AbstractFloat = 1.0E-10, odertol::AbstractFloat = 1.0E-4) where {NS,IntT<:Integer,RealT<:AbstractFloat}

    𝔛 = SparseStateSpace(model.stoich_matrix, initial_distribution.states)
    sink_count = get_sink_count(𝔛)
    p0 = initial_distribution.values

    A = FspMatrixSparse(𝔛, model.propensities, θ = model.parameters)
    u0 = [p0; zeros(sink_count)]
    function odefun!(du, u, θ, t)
        matvec!(du, t, A, u)
        nothing
    end
    fspprob = ODEProblem(odefun!, u0, tspan, p = model.parameters)
    solutions = DE.solve(fspprob, fspalgorithm.ode_method, atol = odeatol, rtol = odertol, saveat = saveat)

    output = SparseFspOutput{NS,IntT,RealT}(
        t = Vector{RealT}(),
        p = Vector{SparseMultIdxVector{NS,IntT,RealT}}(),
        sinks = Vector{Vector{RealT}}()
    )
    for (t, u) in zip(solutions.t, solutions.u)
        push!(output.t, t)
        push!(output.p, SparseMultIdxVector(𝔛.states, u[1:end-sink_count]))
        push!(output.sinks, u[end-sink_count+1:end])
    end

    return output
end

Base.@kwdef mutable struct AdaptiveSparseFsp <: TransientCmeAlgorithm
    ode_method::Union{Nothing,AbstractODEAlgorithm}
    space_adapter::SparseSpaceAdapter
end

function solve(model::CmeModel,
    initial_distribution::SparseMultIdxVector{NS,IntT,RealT},
    tspan::Tuple{AbstractFloat,AbstractFloat},
    fspalgorithm::AdaptiveSparseFsp; saveat = [], fsptol::AbstractFloat = 1.0E-6,
    odeatol::AbstractFloat = 1.0E-10,
    odertol::AbstractFloat = 1.0E-4) where {NS,IntT<:Integer,RealT<:AbstractFloat}

    tstart = min(tspan...)
    tend = max(tspan...)
    adapter = fspalgorithm.space_adapter

    p0 = deepcopy(initial_distribution.values)
    statespace = SparseStateSpace(model.stoich_matrix, initial_distribution.states)
    sink_count = get_sink_count(statespace)
    init!(statespace, adapter, p0, tstart, fsptol)

    tnow = tstart
    unow = [p0; zeros(sink_count)]
    A = FspMatrixSparse{RealT}(statespace, model.propensities, θ = model.parameters)

    output = SparseFspOutput{NS,IntT,RealT}(
        t = Vector{RealT}(),
            p = Vector{SparseMultIdxVector{NS,IntT,RealT}}(),
            sinks = Vector{Vector{RealT}}()
        )
    while tnow < tend
        # Set up callback for checking the growth of FSP error over time
        function fsprhs!(du, u, θ, t)
            matvec!(du, t, A, u)
        end
        function affect!(integrator)
            DE.terminate!(integrator)
        end
        function fsp_error_constraint(u, t, integrator)
            sinks = u[end-sink_count+1:end]            
            return sum(sinks) - fsptol * t / tend
        end
        fsp_cb = DE.ContinuousCallback(
            fsp_error_constraint,
            affect!,
            save_positions = (false, false),
            interp_points = 100,
            abstol = eps()        
            )

        fspprob = DE.ODEProblem(fsprhs!, unow, (tnow, tend), p = model.parameters)
        integrator = DE.init(fspprob, fspalgorithm.ode_method, atol = odeatol, rtol = odertol, callback = fsp_cb, saveat = saveat)

        DE.step!(integrator, tend - tnow, true)

        for (t, u) in zip(integrator.sol.t, integrator.sol.u)
            push!(output.t, t)
            push!(output.p, SparseMultIdxVector(statespace.states, u[1:end-sink_count]))
            push!(output.sinks, u[end-sink_count+1:end])
        end

        tnow = integrator.t
        if tnow < tend
            p = integrator.u[1:end-sink_count]
            sinks = integrator.u[end-sink_count+1:end]
            adapt!(statespace, adapter, p, sinks, tnow, tend, fsptol;integrator=integrator)
            A = FspMatrixSparse{RealT}(statespace, model.propensities, θ = model.parameters)              
            if sum(sinks) >= tnow*fsptol/tend 
                sinks .-= eps()
            end            
            unow = [p; sinks]                              
        else
            u = integrator.u
            push!(output.t, tnow)
            push!(output.p, SparseMultIdxVector(statespace.states, u[1:end-sink_count]))
            push!(output.sinks, u[end-sink_count+1:end])
        end
    end
    return output
end

