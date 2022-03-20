import DifferentialEquations as DE
import Sundials
using DifferentialEquations.DiffEqBase: ODEProblem, AbstractODEAlgorithm
import DifferentialEquations.DiffEqBase: solve

export FixedSparseFsp, AdaptiveFspSparse, solve
include("fspoutput.jl")
include("spaceadapters.jl")

function solve(model::CmeModel,
    initial_distribution::FspVectorSparse{NS,IntT,RealT},
    tspan::Union{Vector,Tuple},
    ode_method::Union{Nothing,AbstractODEAlgorithm}; 
    saveat = [], fsptol::AbstractFloat = 1.0E-6, odeatol::AbstractFloat = 1.0E-10, odertol::AbstractFloat = 1.0E-4) where {NS,IntT<:Integer,RealT<:AbstractFloat}

    statespace = StateSpaceSparse(model.stoich_matrix, initial_distribution.states)
    sink_count = get_sink_count(statespace)
    p0 = initial_distribution.values

    A = FspMatrixSparse(statespace, model.propensities, parameters = model.parameters)
    u0 = [p0; zeros(sink_count)]
    function odefun!(du, u, θ, t)
        matvec!(du, t, A, u)
        nothing
    end
    fspprob = ODEProblem(odefun!, u0, tspan, p = model.parameters)
    solutions = DE.solve(fspprob, ode_method, atol = odeatol, rtol = odertol, saveat = saveat)

    output = FspOutputSparse{NS,IntT,RealT}(
        t = Vector{RealT}(),
        p = Vector{FspVectorSparse{NS,IntT,RealT}}(),
        sinks = Vector{Vector{RealT}}()
    )
    for (t, u) in zip(solutions.t, solutions.u)
        push!(output.t, t)
        push!(output.p, FspVectorSparse(statespace.states, u[1:end-sink_count]))
        push!(output.sinks, u[end-sink_count+1:end])
    end

    return output
end

"""
    mutable struct AdaptiveFspSparse <: TransientCmeAlgorithm

Struct to store adaptive Finite State Projection algorithmic options. This type is intended to work with `solve()` methods that output CME solutions as sparse vectors of type `SparseMultIdxVector`. 

# Fields 

    ode_method
Instance of an ODE algorithm from the `DifferentialEquations.jl` package.

    space_adapter 
Method to adapt the FSP state space when truncation error exceeds user-specified tolerance.

# See also 
[`SparseSpaceAdapter`](@ref), [`SparseMultIdxVector`](@ref), [`CmeModel`](@ref)
"""
Base.@kwdef mutable struct AdaptiveFspSparse <: TransientCmeAlgorithm
    ode_method::Union{Nothing,AbstractODEAlgorithm}
    space_adapter::AbstractSpaceAdapterSparse
end

"""
`solve(model::CmeModel,
initial_distribution::SparseMultIdxVector,
tspan,
fspalgorithm::AdaptiveFspSparse; 
saveat = [], fsptol=1.0E-6,
odeatol=1.0E-10,
odertol=1.0E-4,
verbose=false
)`

Numerical integration for the chemical master equation using an adaptive Finite State Projection (FSP) algorithm. This method is adaptive, meaning that states will be added during the course of integration to ensure the approximation error is below the user-specified tolerance. 
In addition, depending on the input `fspalgorithm`, the space adapter may delete states with low probabilities before expanding the state space. This ensures not only accuracy but also efficiency of the CME integration.

# Arguments 

    model::CmeModel
Chemical Master Equation model.

    initial_distribution::SparseMultIdxVector
Initial distribution.

    tspan::Tuple{AbstractFloat, AbstractFloat}
Timespan for the integration. The first element is the starting time, the second element is the end time of the integration.

    fspalgorithm::AdaptiveFspSparse
An instance of `AdaptiveSparseFsp`, storing information about the specific adaptive FSP method.

    saveat::Vector{AbstractFloat} (optional)
Time points to store the FSP solutions at. Default is `[]`, which means all solutions at every timestep will be stored.

    fsptol (optional)
FSP truncation error tolerance. The FSP state space is adapted so that the output solution vectors sum to greater than or equal to `1-fsptol`. Default: 1.0E-6.

    odertol (optional)
Relative error tolerance for the ODE solver of the ODEs resulting from FSP truncation. Default: 1.0E-4.

    odeatol (optional)
Relative error tolerance for the ODE solver of the ODEs resulting from FSP truncation. Default: 1.0E-10.

    verbose (optional)
Whether to output status when updating the state space. Default: false.

# Returns 
An instance of `FspOutputSparse`.

# See also 
[`CmeModel`](@ref), [`AdaptiveFspSparse`](@ref), [`FspOutputSparse`](@ref).
"""
function solve(model::CmeModel,
    initial_distribution::FspVectorSparse{NS,IntT,RealT},
    tspan::Tuple{AbstractFloat,AbstractFloat},
    fspalgorithm::AdaptiveFspSparse; saveat = [], fsptol::AbstractFloat = 1.0E-6,
    odeatol::AbstractFloat = 1.0E-10,
    odertol::AbstractFloat = 1.0E-4,
    verbose::Bool = false) where {NS,IntT<:Integer,RealT<:AbstractFloat}

    tstart = min(tspan...)
    tend = max(tspan...)
    adapter = fspalgorithm.space_adapter

    p0 = deepcopy(initial_distribution.values)
    statespace = StateSpaceSparse(model.stoich_matrix, initial_distribution.states)
    sink_count = get_sink_count(statespace)
    init!(statespace, adapter, p0, tstart, fsptol)

    tnow = tstart
    unow = [p0; zeros(sink_count)]
    A = FspMatrixSparse{RealT}(statespace, model.propensities, parameters = model.parameters)

    output = FspOutputSparse{NS,IntT,RealT}(
        t = Vector{RealT}(),
            p = Vector{FspVectorSparse{NS,IntT,RealT}}(),
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
        # Error constraint for the intermediate solutions. The intermediate sinks must not grow beyond a linear function of time that reaches `fsptol` at the end time
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
            push!(output.p, FspVectorSparse(statespace.states, u[1:end-sink_count]))
            push!(output.sinks, u[end-sink_count+1:end])
        end

        tnow = integrator.t
        # If the ODE integrator exits before reaching `tend`, it means the intermediate solution with current state space has reached error tolerance threshold
        if tnow < tend             
            p = integrator.u[1:end-sink_count]
            sinks = integrator.u[end-sink_count+1:end]
            # Call adapter to modify the state space and shrink/expand the solution vector as appropriate
            adapt!(statespace, adapter, p, sinks, tnow, tend, fsptol;integrator=integrator)
            A = FspMatrixSparse{RealT}(statespace, model.propensities, parameters = get_parameters(model))              

            # This numerical trick ensures that DifferentialEquations.jl's event detection does not miss the event the updated state space is insufficient
            if sum(sinks) >= tnow*fsptol/tend 
                sinks .-= eps()
            end

            # Extended solution vector to start the ODE integrator again            
            unow = [p; sinks]   
            
            if verbose 
                println("t = $(round(tnow, digits=2)). Update state space. New size: $(get_state_count(statespace)).")
            end
        else # Otherwise, add the final slice to the output
            u = integrator.u
            push!(output.t, tnow)
            push!(output.p, FspVectorSparse(statespace.states, u[1:end-sink_count]))
            push!(output.sinks, u[end-sink_count+1:end])
        end
    end
    return output
end
