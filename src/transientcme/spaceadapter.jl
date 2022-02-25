export SparseSpaceAdapter, SpaceAdapter, adapt!

abstract type SpaceAdapter end
abstract type SparseSpaceAdapter <: SpaceAdapter end 

struct RStepAdapter <: SparseSpaceAdapter
    max_step_count::Int 
end

struct SelectiveRStepAdapter <: SparseSpaceAdapter
    max_step_count::Int
end

function adapt!(statespace::SparseStateSpace, adapter::SparseSpaceAdapter, p::Vector{AbstractFloat}, sinks::Vector{AbstractFloat}, t::AbstractFloat, fsptol::AbstractFloat)
    nothing
end
