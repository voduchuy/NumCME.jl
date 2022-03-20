import LinearAlgebra: mul!, axpy!
import SparseArrays as SArrays

export FspMatrixSparse, get_parameters, get_states, get_rowcount, get_colcount, get_propensities

"""
The transition rate matrix of the finite Markov chain approximation to the CME (including sink states). This type is based on the sparse representation of the truncated state space.
"""
Base.@kwdef mutable struct FspMatrixSparse{NS,NR,IntT<:Integer,RealT<:AbstractFloat} <: FspMatrix
    parameters::Vector{Any}
    states::Vector{MVector{NS,IntT}}

    rowcount::IntT
    colcount::IntT

    propensities::Vector{<:Propensity}

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
get_time(A::FspMatrixSparse) = A.t_cache
get_jointtv_propensity_ids(A::FspMatrixSparse) = A.jointtv_propensity_ids
get_separabletv_propensity_ids(A::FspMatrixSparse) = A.separabletv_propensity_ids
get_timeinvariant_propensity_ids(A::FspMatrixSparse) = A.timeinvariant_propensity_ids

"""
FspMatrixSparse{RealT}(space::AbstractStateSpaceSparse{NS,NR,IntT}, propensity_functions::Vector{<:Propensity}; θ = []) 

Return an instance of `FspMatrixSparse` type with nonzero entries of type `RealT`, based on a sparse state space of type `AbstractStateSpaceSparse` and a vector `propensity_functions` of propensity functions. 

# See also 
`Propensity`
"""
function FspMatrixSparse{RealT}(space::AbstractStateSpaceSparse{NS,NR,IntT,SizeT}, propensity_functions::Vector{<:Propensity}; parameters = []) where {NS,NR,IntT<:Integer,RealT<:AbstractFloat,SizeT<:Integer}
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
        ridxs, cidxs, vs = _generate_sparsematrix_entries(space, propensity_functions[reactionidx].f, reactionidx, parameters, valtype=RealT)
        rowindices = vcat(rowindices, ridxs)
        colindices = vcat(colindices, cidxs)
        vals = vcat(vals, vs)
    end
    timeinvariant_matrix = !isempty(timeinvariant_propensity_ids) ? sparse(rowindices, colindices, vals, rowcount, rowcount) : nothing

    separabletv_factormatrices = Vector{SparseMatrixCSC}()
    for reactionidx in separabletv_propensity_ids
        rowindices, colindices, vals = _generate_sparsematrix_entries(space, propensity_functions[reactionidx].statefactor, reactionidx, parameters, valtype=RealT)
        push!(
            separabletv_factormatrices,
            sparse(rowindices, colindices, vals, rowcount, rowcount)
        )
    end

    jointtv_matrices = Vector{SparseMatrixCSC}()
    for reactionidx in jointtv_propensity_ids
        rowindices, colindices, vals = _generate_sparsematrix_entries(space, nothing, reactionidx, parameters, valtype=RealT)
        push!(
            jointtv_matrices,
            sparse(rowindices, colindices, vals, rowcount, rowcount)
        )
    end

    return FspMatrixSparse{NS,NR,IntT,RealT}(
        parameters = parameters,
        t_cache = -Inf,
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
FspMatrixSparse(space::AbstractStateSpaceSparse, propensity_functions::Vector; θ = []) 

Return an instance of `FspMatrixSparse` type with double-precision(i.e.`Float64`) nonzero entries, based on a sparse state space of type `AbstractStateSpaceSparse` and a vector `propensity_functions` of propensity functions. 
"""
function FspMatrixSparse(space::AbstractStateSpaceSparse{NS,NR,IntT,SizeT}, propensity_functions::Vector{<:Propensity}; parameters = []) where {NS,NR,IntT<:Integer,SizeT<:Integer}
    return FspMatrixSparse{Float64}(space, propensity_functions; parameters)
end


function _generate_sparsematrix_entries(space::AbstractStateSpaceSparse{NS,NR,IntT,SizeT}, statefactor::Union{Nothing,Function}, reactionidx::Integer, θ::Vector = []; valtype::Type{RealT}=Float64) where {NS,NR,IntT<:Integer,RealT<:AbstractFloat,SizeT<:Integer}    
    state_count = get_state_count(space)
    rowindices = zeros(IntT, 2 * state_count)
    colindices = zeros(IntT, 2 * state_count)
    vals = zeros(valtype, 2 * state_count)

    # Not optimal, should evaluate propensity function only once per CME state
    @simd for idx = 1:state_count
        # Populate the diagonal entries
        dval = (statefactor === nothing) ? 0.0 : statefactor(space.states[idx], θ)
        @inbounds rowindices[idx] = idx
        @inbounds colindices[idx] = idx
        @inbounds vals[idx] = -1.0 * dval

        # Populate the row corresponding to a reachable sink state (if any)
        sink_idx = space.sink_connectivity[idx][reactionidx]
        if sink_idx ≠ 0
            @inbounds colindices[state_count+idx] = idx
            @inbounds rowindices[state_count+idx] = state_count + sink_idx
            @inbounds vals[state_count+idx] = dval
        end
    end
    @simd for idx = 1:state_count
        # If nonzero, this will be the column index of the nonzero entry on the idx-th row
        cidx = space.state_connectivity[idx][reactionidx]
        if cidx ≠ 0
            @inbounds rowindices[cidx+state_count] = idx
            @inbounds colindices[cidx+state_count] = cidx
            @inbounds vals[cidx+state_count] = -1.0 * vals[cidx]            
        end
    end
    return rowindices[colindices.≠0], colindices[colindices.≠0], vals[colindices.≠0]
end

function _update_sparsematrix!(matrix::SparseMatrixCSC, states::Vector, propensity::Any, t::AbstractFloat, θ::Vector = [])
    m = matrix.m
    n = length(states)    
    rowvals = SArrays.getrowval(matrix)
    nzvals = SArrays.getnzval(matrix)
    @simd for icol in 1:n
        val = propensity.f(t, states[icol], θ)
        for i in nzrange(matrix, icol)
            @inbounds nzvals[i] = (rowvals[i] == icol) ? -1.0 * val : val
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
"""
matvec!(out::Vector{RealT}, t::AbstractFloat, A::FspMatrixSparse, v::Vector{RealT})

Perform matrix-vector multiplication for the FSP matrix `A` at time `t` and vector `v`. The result is written to the vector `out`.
"""
function matvec!(out, t, A::FspMatrixSparse, v)
    if A.timeinvariant_matrix isa Nothing
        out .= 0.0
    else        
        # out .= A.timeinvariant_matrix*v
        SArrays.mul!(out, A.timeinvariant_matrix, v)
    end
    θ = A.parameters
    for (i, r) in enumerate(A.separabletv_propensity_ids)
        # out .+= A.propensities[r].tfactor(t, θ...)*A.separabletv_factormatrices[i]*v
        SArrays.mul!(out, A.separabletv_factormatrices[i], v, A.propensities[r].tfactor(t, θ), 1)        
    end    
    needupdate = false
    if t ≠ A.t_cache
        needupdate = true
        A.t_cache = t
    end    
    for (i, r) in enumerate(A.jointtv_propensity_ids)
        (needupdate) && _update_sparsematrix!(A.jointtv_matrices[i], A.states, A.propensities[r],
            t, θ)
        # out .+= A.jointtv_matrices[i]*v 
        SArrays.mul!(out, A.jointtv_matrices[i], v, 1.0, 1.0)                        
    end
    return nothing
end

export matvecadd!

"""
    matvecadd!(out, t, A::FspMatrixSparse, v)

Perform the in-place update `out = out + A(t)*v`.
"""
function matvecadd!(out, t, A::FspMatrixSparse, v)
    if !(A.timeinvariant_matrix isa Nothing    )
        SArrays.mul!(out, A.timeinvariant_matrix, v, 1, 1)
    end
    θ = A.parameters
    for (i, r) in enumerate(A.separabletv_propensity_ids)
        # out .+= A.propensities[r].tfactor(t, θ...)*A.separabletv_factormatrices[i]*v
        SArrays.mul!(out, A.separabletv_factormatrices[i], v, A.propensities[r].tfactor(t, θ), 1)        
    end    
    needupdate = false
    if t ≠ A.t_cache
        needupdate = true
        A.t_cache = t
    end    
    for (i, r) in enumerate(A.jointtv_propensity_ids)
        (needupdate) && _update_sparsematrix!(A.jointtv_matrices[i], A.states, A.propensities[r],
            t, θ)
        # out .+= A.jointtv_matrices[i]*v 
        SArrays.mul!(out, A.jointtv_matrices[i], v, 1.0, 1.0)                        
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

import Base: *
function *(A_at_t::Tuple{AbstractFloat,FspMatrixSparse}, v::Vector{RealT}) where {RealT<:AbstractFloat}
    t = A_at_t[1]
    A = A_at_t[2]
    return matvec(t, A, v)
end

function *(A::FspMatrixSparse, v::Vector{RealT}) where {RealT<:AbstractFloat}
    return matvec(0.0, A, v)
end

