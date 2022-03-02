export CmeModel, CmeModelWithSensitivity, get_stoich_matrix, get_propensities, gradient_sparsity_patterns, propensity_gradients

include("propensity.jl")

struct CmeModel{IntT<:Integer,PT<:Propensity}
    stoich_matrix::Matrix{IntT}
    propensities::Vector{PT}
end
# Accessors 
get_stoich_matrix(model::CmeModel) = model.stoich_matrix
get_propensities(model::CmeModel) = model.propensities

struct CmeModelWithSensitivity{IntT<:Integer,PT<:Propensity,PGradT<:PropensityGradient}
    cmemodel::CmeModel
    gradient_sparsity_patterns::SparseMatrixCSC{IntT}
    propensity_gradients::Vector{PGradT}
end
gradient_sparsity_patterns(model::CmeModelWithSensitivity) = model.gradient_sparsity_patterns
propensity_gradients(model::CmeModelWithSensitivity) = model.propensity_gradients
@forward((CmeModelWithSensitivity, :cmemodel), CmeModel)
