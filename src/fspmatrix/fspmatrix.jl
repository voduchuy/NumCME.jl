using SparseArrays
import Base: *, size
import LinearAlgebra: mul!, axpy!

export FspMatrix, FspSparseMatrix, matvec, matvec!

abstract type FspMatrix end


Base.@kwdef mutable struct FspSparseMatrix{NS,NR,IntT<:Integer,RealT<:AbstractFloat} <: FspMatrix
    space::AbstractSparseStateSpace{NS,NR,IntT}
    rowcount::IntT
    colcount::IntT

    propensity_functions::Vector{TimeSeparablePropensity}
    tvarying_reactions::Vector{IntT}

    tinvariant_matrix::Any
    tvarying_matrices::Vector{Any}

    workvec::Vector{RealT}
end

function FspSparseMatrix{RealT}(space::AbstractSparseStateSpace{NS,NR,IntT}, propensity_functions::Vector{TimeSeparablePropensity}) where {NS,NR,IntT<:Integer,RealT<:AbstractFloat}
    state_count = get_state_count(space)
    sink_count = get_sink_count(space)

    tvarying_reactions = Vector{IntT}()
    tinvariant_reactions = Vector{IntT}()
    for (ir, propensity) in enumerate(propensity_functions)
        push!((propensity.tfactor ≠ nothing) ? tvarying_reactions : tinvariant_reactions, ir)
    end

    n = state_count + sink_count
    tvarying_matrices = Vector()
    for reactionidx in tvarying_reactions
        rowindices, colindices, vals = _generate_sparsematrix_entries(space, propensity_functions[reactionidx].statefactor, reactionidx)
        push!(
            tvarying_matrices,
            sparse(rowindices, colindices, vals, n, n)
        )
    end

    rowindices = Vector{IntT}()
    colindices = Vector{IntT}()
    vals = Vector{RealT}()
    for reactionidx in tinvariant_reactions
        ridxs, cidxs, vs = _generate_sparsematrix_entries(space, propensity_functions[reactionidx].statefactor, reactionidx)
        rowindices = vcat(rowindices, ridxs)
        colindices = vcat(colindices, cidxs)
        vals = vcat(vals, vs)
    end
    tinvariant_matrix = !isempty(tinvariant_reactions) ? sparse(rowindices, colindices, vals, n, n) : nothing

    return FspSparseMatrix{NS,NR,IntT,RealT}(
        space = space,
        rowcount = n,
        colcount = n,
        propensity_functions = propensity_functions,
        tvarying_reactions = tvarying_reactions,
        tinvariant_matrix = tinvariant_matrix,
        tvarying_matrices = tvarying_matrices,
        workvec = zeros(RealT, n)
    )
end

function FspSparseMatrix(space::AbstractSparseStateSpace{NS,NR,IntT}, propensity_functions::Vector{TimeSeparablePropensity}) where {NS,NR,IntT<:Integer}
    return FspSparseMatrix{Float64}(space, propensity_functions)
end


function _generate_sparsematrix_entries(space::AbstractSparseStateSpace{NS,NR,IntT}, statefactor::Function, reactionidx::Integer) where {NS,NR,IntT<:Integer}
    state_count = get_state_count(space)

    rowindices = Vector{IntT}()
    colindices = Vector{IntT}()
    vals = Vector{Float64}()

    sizehint!(rowindices, 2 * state_count)
    sizehint!(colindices, 2 * state_count)
    sizehint!(vals, 2 * state_count)

    for idx = 1:state_count

        # If nonzero, this will be the column index of the nonzero entry on the idx-th row
        cidx = space.state_connectivity[idx][reactionidx]
        if cidx ≠ 0
            push!(rowindices, idx)
            push!(colindices, cidx)
            push!(vals, statefactor(space.states[cidx]...))
        end

        # Populate the diagonal entries
        dval = statefactor(space.states[idx]...)
        push!(rowindices, idx)
        push!(colindices, idx)
        push!(vals, -1.0 * dval)

        # Populate the row corresponding to a reachable sink state (if any)
        sink_idx = space.sink_connectivity[idx][reactionidx]
        if sink_idx ≠ 0
            push!(rowindices, state_count + sink_idx)
            push!(colindices, idx)
            push!(vals, dval)
        end
    end
    return rowindices[colindices.≠0], colindices[colindices.≠0], vals[colindices.≠0]
end

# Functions to query the dimensions of the matrix 
function size(A::FspSparseMatrix)
    return (A.rowcount, A.colcount)
end

function size(A::FspSparseMatrix, dim::Integer)
    if !(1 ≤ dim ≤ 2)
        throw(ArgumentError("Second argument must be either 1 or 2."))
    end
    return (dim == 1) ? A.rowcount : A.colcount
end

# Matrix-vector multiplications
function matvec(t::AbstractFloat, A::FspSparseMatrix, v::Vector{RealT}) where {RealT <: AbstractFloat}
    w = !(A.tinvariant_matrix isa Nothing) ? A.tinvariant_matrix * v : zeros(RealT, length(v))
    for (i, r) in enumerate(A.tvarying_reactions)
        w += A.propensity_functions[r].tfactor(t)*A.tvarying_matrices[i]*v
        # The following should have saved time and allocations but somehow it INCREASES both!!!
        # mul!(w, A.tvarying_matrices[i], v, A.propensity_functions[r].tfactor(t), 1.0)        
    end
    return w
end

function matvec!(t::AbstractFloat, A::FspSparseMatrix, v::Vector{RealT}, out::Vector{RealT}) where {RealT <: AbstractFloat}    
    if A.tinvariant_matrix isa Nothing 
        out .= 0.0
    else         
        mul!(out, A.tinvariant_matrix, v)
    end

    for (i, r) in enumerate(A.tvarying_reactions)
        out += A.propensity_functions[r].tfactor(t)*A.tvarying_matrices[i]*v
        # The following should have saved time and allocations but somehow it INCREASES both!!!
        # mul!(out, A.tvarying_matrices[i], v, A.propensity_functions[r].tfactor(t), 1.0)
        # Same problem with this
        # axpy!(A.propensity_functions[r].tfactor(t), A.workvec, out)
    end
    return nothing
end

function (A::FspSparseMatrix)(t::AbstractFloat)
    return (t, A)
end

function *(A_at_t::Tuple{AbstractFloat,FspSparseMatrix}, v::Vector{RealT}) where {RealT <: AbstractFloat}
    t = A_at_t[1]
    A = A_at_t[2]
    return matvec(t, A, v)
end

function *(A::FspSparseMatrix, v::Vector{RealT}) where {RealT <: AbstractFloat}
    return matvec(0.0, A, v)
end
