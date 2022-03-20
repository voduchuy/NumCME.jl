export CmeModel, CmeModelWithSensitivity, get_stoich_matrix, get_propensities, gradient_sparsity_patterns, propensity_gradients, get_species_count, get_reaction_count, get_parameter_count, get_parameters
export get_propensity_gradients, get_gradient_sparsity_patterns

include("propensity.jl")
include("propensitygrad.jl")

"""
Stochastic reaction network model.

## Fields 
    stoich_matrix::Matrix{IntT}
(Net) stoichiometry matrix. Column `i` represents the net change to molecular counts due to reaction `i`.

    propensities::Vector{PT}
Propensity functions. 

    parameters::AbstractVector
Parameters.

## Examples

The following code constructs a CmeModel instance for the two-state telegraphic gene expression model.
    
```julia
𝕊 = [[-1, 1, 0] [1, -1, 0] [0, 0, 1] [0, 0, -1]] # Net stoichiometry matrix for 3 species `inactive_gene`, `activated_gene`, `mRNA`
a1 = propensity() do x, p
    p[1] * x[1]
end
a2 = propensity() do x, p
    p[2] * x[2]
end
a3 = propensity() do x, p
    p[3] * x[2]
end
a4 = propensity() do x, p
    p[4] * x[3]
end
k₀₁ = 0.05
k₁₀ = 0.1
λ = 5.0
γ = 0.5
θ = [k₀₁, k₁₀, λ, γ]
model = CmeModel(𝕊, [a1,a2,a3,a4], θ)
```

## See also 
[`CmeModelWithSensitivity`](@ref), [`Propensity`](@ref), [`StandardTimeInvariantPropensity`](@ref), [`JointTimeVaryingPropensity`](@ref), [`SeparableTimeVaryingPropensity`](@ref).
"""
struct CmeModel{IntT<:Integer}
    stoich_matrix::Matrix{IntT}
    propensities::Vector{<:Propensity}
    parameters::AbstractVector
end
# Accessors 
"""
Get the vector of model parameters.
"""
get_parameters(model::CmeModel) = model.parameters
"""
Get number of species.
"""
get_species_count(model::CmeModel) = size(model.stoich_matrix, 1)
"""
Get number of reactions.
"""
get_reaction_count(model::CmeModel) = size(model.stoich_matrix, 2)
"""
Get number of parameters.
"""
get_parameter_count(model::CmeModel) = length(model.parameters)
"""
Get net stoichiometry matrix.
"""
get_stoich_matrix(model::CmeModel) = model.stoich_matrix
"""
Get vector of propensities.
"""
get_propensities(model::CmeModel) = model.propensities

# Pretty printing for CmeModel
using Printf 
function Base.show(io::IO, model::CmeModel)
    println(io, "Stochastic reaction network with $(get_species_count(model)) species, $(get_reaction_count(model)) reactions and $(get_parameter_count(model)) parameters.")
    println(io, "Net stoichiometry matrix")
    Base.print_matrix(io, model.stoich_matrix)
    println(io)
    nothing 
end

"""
    CmeModelWithSensitivity

Stochastic reaction network model with information about gradient of the propensity functions with respect to parameters. 

## See also 
[`PropensityGradient`](@ref)
"""
mutable struct CmeModelWithSensitivity{IntT<:Integer}
    cmemodel::CmeModel{IntT}
    gradient_sparsity_patterns::SparseMatrixCSC
    propensity_gradients::Vector{<:PropensityGradient}
end
# Accessors 
"""
Get parameter vectors of a Chemical Master Equation model instance.
"""
get_parameters(model::CmeModelWithSensitivity) = get_parameters(model.cmemodel)
"""
Get the number of species (i.e., the length of the CME state vector).
"""
get_species_count(model::CmeModelWithSensitivity) = get_species_count(model.cmemodel)
"""
Get number of reactions in the model.
"""
get_reaction_count(model::CmeModelWithSensitivity) = get_reaction_count(model.cmemodel)
"""
Get number of parrameters.
"""
get_parameter_count(model::CmeModelWithSensitivity) = get_parameter_count(model.cmemodel)
"""
Get the net stoichiometry matrix.
"""
get_stoich_matrix(model::CmeModelWithSensitivity) = get_stoich_matrix(model.cmemodel)
"""
Get the vector of propensities.
"""
get_propensities(model::CmeModelWithSensitivity) = get_propensities(model.cmemodel)
"""
Get the sparse matrix representing dependency pattern of propensity functions of model parameters.
"""
get_gradient_sparsity_patterns(model::CmeModelWithSensitivity) = model.gradient_sparsity_patterns
"""
Get the propensity gradients.
"""
get_propensity_gradients(model::CmeModelWithSensitivity) = model.propensity_gradients


"""
    CmeModelWithSensitivity(model::CmeModel)

Construct a CME model with sensitivity from a CME model (without sensitivity information) using `ForwardDiff.jl` automatic differentiation. 
"""
function CmeModelWithSensitivity(model::CmeModel)
    species_count = get_species_count(model)
    parameter_count = get_parameter_count(model)
    gradient_sparsity_patterns = propensitygrad_sparsity_pattern(species_count, parameter_count, get_propensities(model), get_parameters(model))
    propensity_gradients = Vector{PropensityGradient}()
    for propensity in get_propensities(model)
        push!(propensity_gradients, propensity_forwarddiff(propensity, parameter_count))
    end
    return CmeModelWithSensitivity(
        model,
        gradient_sparsity_patterns,
        propensity_gradients
    )
end

"""
    CmeModelWithSensitivity(stoich_matrix, propensities::Vector{PT}, parameters::AbstractVector) where {PT <: Propensity}

Construct a CME model with sensitivity information using a net stoichiometry matrix `stoich_matrix` and propensity functions `propensities` and parameters `parameters`. The partial derivatives of the propensity functions with respect to parameters are automatically derived using automatic differentiation.
"""
CmeModelWithSensitivity(stoich_matrix, propensities::Vector{<:Propensity}, parameters::AbstractVector) where {PT <: Propensity} = CmeModelWithSensitivity(CmeModel(stoich_matrix, propensities, parameters))

# Pretty printing for CmeModelWithSensitivity
using Printf 
function Base.show(io::IO, model::CmeModelWithSensitivity)
    Base.show(model.cmemodel)
    println(io, "Propensity gradient sparsity pattern (reaction × parameter):")
    Base.print_matrix(io, get_gradient_sparsity_patterns(model))
    nothing 
end

include("catalyst_interface.jl")