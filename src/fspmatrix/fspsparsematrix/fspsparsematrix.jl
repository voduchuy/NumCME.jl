using EtFsp
using SparseArrays
using StaticArrays
import LinearAlgebra: mul!, axpy!

export FspMatrixSparse, get_parameters, get_states, get_rowcount, get_colcount, get_propensities 

"""
The transition rate matrix of the finite Markov chain approximation to the CME (including sink states). This type is based on the sparse representation of the truncated state space.
"""
Base.@kwdef mutable struct FspMatrixSparse{NS,NR, PT<: Propensity, IntT<:Integer,RealT<:AbstractFloat} <: FspMatrix
    parameters::Vector{Any}     
    states::Vector{MVector{NS, IntT}}

    rowcount::IntT
    colcount::IntT
    
    propensities::Vector{PT}

    timeinvariant_propensity_ids::Vector{IntT}
    timeinvariant_matrix::SparseMatrixCSC
    
    t_cache::RealT 
    jointtv_propensity_ids::Vector{Int}
    jointtv_matrices::Vector{SparseMatrixCSC}
    
    separabletv_propensity_ids::Vector{IntT}
    separabletv_factormatrices::Vector{SparseMatrixCSC}    
end
# Accessors 
get_parameters(A::FspMatrixSparse) = A.parameters 
get_states(A::FspMatrixSparse) = A.states 
get_rowcount(A::FspMatrixSparse) = A.rowcount 
get_colcount(A::FspMatrixSparse) = A.colcount 
get_propensities(A::FspMatrixSparse) = A.propensities 

"""
FspMatrixSparse{RealT}(space::AbstractSparseStateSpace{NS,NR,IntT}, propensity_functions::Vector{PT}; θ = []) 

Return an instance of `FspMatrixSparse` type with nonzero entries of type `RealT`, based on a sparse state space of type `AbstractSparseStateSpace` and a vector `propensity_functions` of propensity functions. 
"""
function FspMatrixSparse{RealT}(space::AbstractSparseStateSpace{NS,NR,IntT}, propensity_functions::Vector{PT}; θ = []) where {NS,NR,PT<:Propensity, IntT<:Integer,RealT<:AbstractFloat}
    state_count = get_state_count(space)
    sink_count = get_sink_count(space)
    rowcount = state_count + sink_count
    colcount = rowcount 

    timeinvariant_propensity_ids = Vector{IntT}()
    jointtv_propensity_ids = Vector{IntT}()
    separabletv_propensity_ids = Vector{IntT}()
    
    for (ir, propensity) in enumerate(propensity_functions)
        if !istimevarying(propensity)
            push!(timeinvariant_propensity_ids, ir)
        else 
            push!(istimeseparable(propensity) ? separabletv_propensity_ids : jointtv_propensity_ids, ir)
        end
    end

    rowindices = Vector{IntT}()
    colindices = Vector{IntT}()
    vals = Vector{RealT}()
    for reactionidx in timeinvariant_propensity_ids
        ridxs, cidxs, vs = _generate_sparsematrix_entries(space, propensity_functions[reactionidx].f, reactionidx, θ = θ)
        rowindices = vcat(rowindices, ridxs)
        colindices = vcat(colindices, cidxs)
        vals = vcat(vals, vs)
    end
    timeinvariant_matrix = !isempty(timeinvariant_propensity_ids) ? sparse(rowindices, colindices, vals, rowcount, rowcount) : nothing

    separabletv_factormatrices = Vector{SparseMatrixCSC}()
    for reactionidx in separabletv_propensity_ids 
        rowindices, colindices, vals = _generate_sparsematrix_entries(space, propensity_functions[reactionidx].statefactor, reactionidx, θ = θ)
        push!(
            separabletv_factormatrices,
            sparse(rowindices, colindices, vals, rowcount, rowcount)
        )
    end

    jointtv_matrices = Vector{SparseMatrixCSC}()
    for reactionidx in jointtv_propensity_ids
        rowindices, colindices, vals = _generate_sparsematrix_entries(space, nothing, reactionidx, θ = θ)
        push!(
            jointtv_matrices,
            sparse(rowindices, colindices, vals, rowcount, rowcount)
        )
    end

    return FspMatrixSparse{NS,NR,PT,IntT,RealT}(
        parameters = θ,
        t_cache = 0.0,
        states = deepcopy(space.states),
        rowcount = rowcount,
        colcount = colcount,
        propensities = propensity_functions,        
        timeinvariant_propensity_ids = timeinvariant_propensity_ids,
        timeinvariant_matrix = timeinvariant_matrix,
        jointtv_propensity_ids = jointtv_propensity_ids,
        jointtv_matrices = jointtv_matrices,
        separabletv_propensity_ids = separabletv_propensity_ids,
        separabletv_factormatrices = separabletv_factormatrices        
    )
end

"""
FspMatrixSparse(space::AbstractSparseStateSpace, propensity_functions::Vector; θ = []) 

Return an instance of `FspMatrixSparse` type with double-precision(i.e.`Float64`) nonzero entries, based on a sparse state space of type `AbstractSparseStateSpace` and a vector `propensity_functions` of propensity functions. 
"""
function FspMatrixSparse(space::AbstractSparseStateSpace{NS,NR,IntT}, propensity_functions::Vector{PT}; θ  =[]) where {NS,NR,PT<:Propensity,IntT<:Integer}
    return FspMatrixSparse{Float64}(space, propensity_functions; θ)
end


function _generate_sparsematrix_entries(space::AbstractSparseStateSpace{NS,NR,IntT}, statefactor::Union{Nothing, Function}, reactionidx::Integer; θ::Vector = []) where {NS,NR,IntT<:Integer}
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
            push!(vals, (statefactor === nothing) ? 0.0 : statefactor(space.states[cidx]..., θ...))
        end

        # Populate the diagonal entries
        dval = (statefactor === nothing) ? 0.0 : statefactor(space.states[idx]..., θ...)
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

function _update_sparsematrix!(matrix::SparseMatrixCSC, states::Vector, propensity::JointTimeVaryingPropensity, t::AbstractFloat, θ::Vector = []) 
    m = matrix.m 
    n = length(states) 
    colptr = matrix.colptr 
    rowval = matrix.rowval
    nzval = matrix.nzval 
    for icol in 1:n
        val = propensity.f(t, states[icol]..., θ...)
        for i in colptr[icol]:(colptr[icol+1]-1)
           irow = rowval[i]
           nzval[i] = (irow == icol) ? -1.0*val : val
        end
    end
    nothing
end

# Functions to query the dimensions of the matrix 
import Base: size 

"""
Return the size of the FSP matrix.
"""
function size(A::FspMatrixSparse)
    return (A.rowcount, A.colcount)
end

"""
Return number of rows (if `dim=1`) or columns (if `dim=2`) of the FSP matrix.
"""
function size(A::FspMatrixSparse, dim::Integer)
    if !(1 ≤ dim ≤ 2)
        throw(ArgumentError("Second argument must be either 1 or 2."))
    end
    return (dim == 1) ? A.rowcount : A.colcount
end

# Matrix-vector multiplications
import Base: *

"""
matvec!(out::Vector{RealT}, t::AbstractFloat, A::FspMatrixSparse, v::Vector{RealT})

Perform matrix-vector multiplication for the FSP matrix `A` at time `t` and vector `v`. The result is written to the vector `out`.
"""
function matvec!(out, t, A::FspMatrixSparse, v) 
    if A.timeinvariant_matrix isa Nothing
        out .= 0.0
    else
        mul!(out, A.timeinvariant_matrix, v)
    end
    θ = A.parameters    
    for (i, r) in enumerate(A.separabletv_propensity_ids)
        out += A.propensities[r].tfactor(t, θ...) * A.separabletv_factormatrices[i] * v
        # Using mul!() should have saved time and allocations but somehow it INCREASES both!!!
    end
    (isempty(A.jointtv_propensity_ids)) && return nothing
    needupdate = false 
    if t ≠ A.t_cache
        needupdate = true 
        A.t_cache = t 
    end 
    for (i, r) in enumerate(A.jointtv_propensity_ids)
        (needupdate) && _update_sparsematrix!(A.jointtv_matrices[i], A.states, A.propensities[r],
        t, θ)
        out += A.jointtv_matrices[i] * v
    end
    return nothing
end

"""
w = matvec(out::Vector{RealT}, t::AbstractFloat, A::FspMatrixSparse, v::Vector{RealT})

Perform matrix-vector multiplication for the FSP matrix `A` at time `t` and vector `v` and return `w` as the result.
"""
function matvec(t, A::FspMatrixSparse, v) 
    w = similar(v)
    matvec!(w, t, A, v)
    return w
end

function (A::FspMatrixSparse)(t::AbstractFloat)
    return (t, A)
end

function *(A_at_t::Tuple{AbstractFloat,FspMatrixSparse}, v::Vector{RealT}) where {RealT<:AbstractFloat}
    t = A_at_t[1]
    A = A_at_t[2]
    return matvec(t, A, v)
end

function *(A::FspMatrixSparse, v::Vector{RealT}) where {RealT<:AbstractFloat}
    return matvec(0.0, A, v)
end

