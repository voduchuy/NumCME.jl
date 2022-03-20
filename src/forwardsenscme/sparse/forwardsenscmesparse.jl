export AdaptiveForwardSensFspSparse, ForwardSensFspInitialConditionSparse
export forwardsens_initial_condition, get_states, get_probability, get_sensitivity

"""
    $(TYPEDEF)

Type to store algorithmic options for solving a Foward Sensitivity CME problem.

# Fields 

$(FIELDS)

# See also 
[`ForwardSensRStepAdapter`](@ref)
"""
Base.@kwdef struct AdaptiveForwardSensFspSparse <: ForwardSensCmeAlgorithm
    """
    Method to adapt state space when the current solution error exceeds tolerance.
    """
    space_adapter::AbstractForwardSensSpaceAdapterSparse
    """
    An ODE solver from `DifferentialEquations.jl` to integrate the FSP-truncated problem.
    """
    ode_method::Union{Nothing,AbstractODEAlgorithm}
end

"""
    $(TYPEDEF)

Initial condition for Foward Sensitivity FSP in sparse format. 

# Fields 

$(FIELDS)

"""
struct ForwardSensFspInitialConditionSparse{NS,IntT<:Integer,RealT<:AbstractFloat}
    """
    List of states with nonzero initial probabilities.
    """
    states::Vector{MVector{NS,IntT}}
    """
    List of nonzero probabilities, with `p[i]` being the probability of state `states[i]`.
    """
    p::Vector{RealT}
    """
    List of partial derivatives, with `S[j][i]` being the `j`-th partial derivative of `p[i]`. 
    """
    S::Vector{Vector{RealT}}
end
get_states(ic::ForwardSensFspInitialConditionSparse) = ic.states
get_probability(ic::ForwardSensFspInitialConditionSparse) = ic.p
get_sensitivity(ic::ForwardSensFspInitialConditionSparse) = ic.S

"""
    $(TYPEDSIGNATURES)

Return a `ForwardSensFspInitialConditionSparse` instance from a list of initial states `states`, initial probabilities `probabilities` and initial sensitivities `sensitivity`.   
"""
function forwardsens_initial_condition(states::Vector{<:AbstractVector{IntT}}, probabilities::Vector{RealT}, sensitivity::Vector{Vector{RealT}}) where {IntT<:Integer, RealT<:AbstractFloat}
    isempty(states) && throw(ArgumentError("Empty state list in input."))
    species_count = length(states[1])    
    
    _states = Vector{MVector{species_count, IntT}}()
    for eachstate in states 
        push!(_states, eachstate)
    end
    _p = copy(probabilities)
    _S = deepcopy(sensitivity)
    return ForwardSensFspInitialConditionSparse{species_count, IntT, RealT}(
        _states, 
        _p, 
        _S 
    )
end

"""
$(TYPEDSIGNATURES)

Foward sensitivity Finite State Projection to compute the time-dependent probability distributions of a stochastic reaction network along with their partial derivatives with respect to model parameters.

# Arguments 
- `model`: Stochastic reaction network model with propensity sensitivities.
- `initial_condition`: Initial condition.
- `tspan`: Time span. 
- `sensfspalgorithm`: Forward Sesitivity FSP method, a combination of space adaptation strategy and the ODE method to integrate the FSP-truncated problem.
- `saveat`: timepoints to save the intermediate solutions. If left empty, all timesteps computed by the ODE solver is saved.
- `fsptol`: FSP error tolerance. The solver dynamically adapts the projection space so that the sum of probability mass in the truncated solution is above `1.0 - fsptol`. Default: 1.0E-6.
- `odeatol`: Absolute tolerance for the ODE solver. Default: 1.0E-10.
- `odertol`: Relative tolerance for the ODE solver. Default: 1.0E-4.
- `verbose`: Whether to print status when adapting state space. Default: false.

# Returns 
A `ForwardSensFspOutputSparse` instance.

# See also 
[`FspVectorSparse`](@ref), [`ForwardSensFspOutputSparse`](@ref), [`AdaptiveForwardSensFspSparse`](@ref).
"""
function solve(model::CmeModelWithSensitivity,
    initial_condition::ForwardSensFspInitialConditionSparse{NS, IntT, RealT},
    tspan::Tuple{AbstractFloat,AbstractFloat},
    sensfspalgorithm::AdaptiveForwardSensFspSparse;
    saveat = [],
    fsptol::AbstractFloat = 1.0E-6,
    odeatol::AbstractFloat = 1.0E-10,
    odertol::AbstractFloat = 1.0E-4,
    verbose::Bool=false) where {NS, IntT<:Integer, RealT<:AbstractFloat}

    tstart = min(tspan...)
    tend = max(tspan...)
    parameter_count = get_parameter_count(model)

    (length(get_sensitivity(initial_condition)) ≠ parameter_count) && throw(ArgumentError("Initial condition does not match CME model. Initial condition must contain `np` sensitivity vectors where `np` is the number of CME model parameters."))


    adapter = sensfspalgorithm.space_adapter
    statespace = StateSpaceSparse(get_stoich_matrix(model), get_states(initial_condition))
    sink_count = get_sink_count(statespace)
    p0 = copy(get_probability(initial_condition))
    S0 = deepcopy(get_sensitivity(initial_condition))
    init!(statespace, adapter, p0, S0, tstart, fsptol)


    output = ForwardSensFspOutputSparse{NS,IntT,RealT}(
        t = Vector{RealT}(),
        p = Vector{FspVectorSparse{NS,IntT,RealT}}(),
        sinks = Vector{Vector{RealT}}(),
        S = Vector{Vector{FspVectorSparse{NS,IntT,RealT}}}(),
        dsinks = Vector{Vector{Vector{RealT}}}()
    )        
    
    tnow = tstart
    n = length(p0)+sink_count
    unow = zeros(RealT, n*(parameter_count+1))
    copyto!(unow, 1, p0, 1, length(p0))
    for i in 1:parameter_count
        copyto!(unow, i*n+1, S0[i], 1, length(S0[i]))
    end
    while tnow < tend
        SA = ForwardSensFspMatrixSparse{RealT}(model, statespace)
        n = get_rowcount(SA.fspmatrix)
        function sensfsprhs!(du, u, θ, t)
            matvec!(du, t, SA, u)
        end
        function affect!(integrator)
            DE.terminate!(integrator)
        end
        # Error constraint for the intermediate solutions. The intermediate sinks must not grow beyond a linear function of time that reaches `fsptol` at the end time
        function fsp_error_constraint(u, t, integrator)
            sinks = u[n-sink_count+1:n]
            return sum(sinks) - fsptol * t / tend
        end

        fsp_cb = DE.ContinuousCallback(
            fsp_error_constraint,
            affect!,
            save_positions = (false, false),
            interp_points = 100,
            abstol = eps()
        )

        fsensfspprob = DE.ODEProblem(sensfsprhs!, unow, (tnow, tend), p = get_parameters(model))
        integrator = DE.init(fsensfspprob, sensfspalgorithm.ode_method, atol = odeatol, rtol = odertol, callback = fsp_cb, saveat = saveat, progress=true)
        DE.step!(integrator, tend - tnow, true)        

        for (t, u) in zip(integrator.sol.t, integrator.sol.u)
            push!(output.t, t)
            push!(output.p, FspVectorSparse(statespace.states, u[1:n-sink_count]))
            push!(output.sinks, u[n-sink_count+1:n])
            push!(output.S, [FspVectorSparse(get_states(statespace), u[ip*n+1:(ip+1)*n-sink_count]) for ip in 1:parameter_count])
            push!(output.dsinks, [u[(ip+1)*n-sink_count+1:(ip+1)*n] for ip in 1:parameter_count])
        end

        tnow = integrator.t
        # If the ODE integrator exits before reaching `tend`, it means the intermediate solution with current state space has reached error tolerance threshold
        if tnow < tend            
            n = get_rowcount(SA.fspmatrix)
            p = integrator.u[1:n-sink_count]
            sinks = integrator.u[n-sink_count+1:n]
            S = [integrator.u[ip*n+1:(ip+1)*n-sink_count] for ip in 1:parameter_count]
            dsinks = [integrator.u[(ip+1)*n-sink_count+1:(ip+1)*n] for ip in 1:parameter_count]
            
            # Call adapter to modify the state space and shrink/expand the solution vector as appropriate
            adapt!(statespace, adapter, p, S, sinks, dsinks, tnow, tend, fsptol, integrator=integrator)

            # This numerical trick ensures that DifferentialEquations.jl's event detection does not miss the event the updated state space is insufficient
            if sum(sinks) >= tnow * fsptol / tend
                sinks .-= eps()
            end

            # Extended solution vector to start the ODE integrator again            
            n = length(p)+sink_count
            unow = zeros(RealT, n*(parameter_count+1))
            copyto!(unow, 1, p, 1, length(p))
            copyto!(unow, n - sink_count+1, sinks, 1, sink_count)
            for ip in 1:parameter_count
                copyto!(unow, ip*n+1, S[ip], 1, length(S[ip]))
                copyto!(unow, (ip+1)*n - sink_count + 1, dsinks[ip], 1, sink_count)
            end
            if verbose 
                println("At t = $(round(tnow, digits=2)): update sate space. New size: $(get_state_count(statespace)) states.")
            end
        else # Otherwise, add the final slice to the output
            u = integrator.u            
            push!(output.t, tnow)
            push!(output.p, FspVectorSparse(statespace.states, u[1:n-sink_count]))
            push!(output.sinks, u[n-sink_count+1:n])
            push!(output.S, [FspVectorSparse(get_states(statespace), u[ip*n+1:ip*n+n-sink_count]) for ip in 1:parameter_count])
            push!(output.dsinks, [u[ip*n+n-sink_count+1:ip*n+n] for ip in 1:parameter_count])
        end
    end
    return output
end