export CmeModel, CmeModelWithSensitivity, get_stoich_matrix, get_propensities, gradient_sparsity_patterns, propensity_gradients, get_species_count, get_reaction_count, get_parameter_count
export get_propensity_gradients, get_gradient_sparsity_patterns

include("propensity.jl")
include("propensitygrad.jl")

struct CmeModel{IntT<:Integer,PT<:Propensity}
    stoich_matrix::Matrix{IntT}
    propensities::Vector{PT}
    parameters::AbstractVector
end
# Accessors 
get_species_count(model::CmeModel) = size(model.stoich_matrix, 1)
get_reaciton_count(model::CmeModel) = size(model.stoich_matrix, 2)
get_parameter_count(model::CmeModel) = length(model.parameters)
get_stoich_matrix(model::CmeModel) = model.stoich_matrix
get_propensities(model::CmeModel) = model.propensities

mutable struct CmeModelWithSensitivity{IntT<:Integer,PT<:Propensity,PGradT<:PropensityGradient}
    cmemodel::CmeModel{IntT, PT}
    gradient_sparsity_patterns::SparseMatrixCSC
    propensity_gradients::Vector{PGradT}
end
# Accessors 
get_gradient_sparsity_patterns(model::CmeModelWithSensitivity) = model.gradient_sparsity_patterns
get_propensity_gradients(model::CmeModelWithSensitivity) = model.propensity_gradients
@forward((CmeModelWithSensitivity, :cmemodel), CmeModel)


# Generate a CME model with sensitivity from a CME model using ForwardDiff automatic differentiation 
function CmeModelWithSensitivity(model::CmeModel)
    species_count = get_species_count(model)
    parameter_count = get_parameter_count(model)
    gradient_sparsity_patterns = propensitygrad_sparsity_pattern(species_count, parameter_count, get_propensities(model))
    propensity_gradients = Vector{PropensityGradient}()
    for propensity in get_propensities(model)
        push!(propensity_gradients, propensity_forwarddiff(propensity, species_count, parameter_count))
    end
    return CmeModelWithSensitivity(
        model,
        gradient_sparsity_patterns,
        propensity_gradients
    )
end
