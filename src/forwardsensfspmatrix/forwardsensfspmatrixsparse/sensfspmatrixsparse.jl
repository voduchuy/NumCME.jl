export AbstractForwardSensFspMatrixSparse, ForwardSensFspMatrixSparse, matvec!

abstract type AbstractForwardSensFspMatrixSparse <: ForwardSensFspMatrix end

"""
The time-varying matrix that defines the right-hand side of the truncated forward sensitivity Chemical Master Equation. 
This type works with `StateSpaceSparse` subtypes and it uses `SparseMatrixCSC` type to internally represent component matrices.
"""
mutable struct ForwardSensFspMatrixSparse{NS,NR,IntT<:Integer,RealT<:AbstractFloat} <: AbstractForwardSensFspMatrixSparse
    fspmatrix::FspMatrixSparse{NS,NR,IntT,RealT}

    propensity_gradients::Vector{<:PropensityGradient}

    timeinvariant_matdiffs::Vector{SparseMatrixCSC}

    separabletv_matdiff_sparsity_pattern::SparseMatrixCSC{Bool}
    separabletv_nzmatdiffs::Vector{SparseMatrixCSC}
    separabletv_mat_map::Vector{IntT}

    jointtv_matdiff_sparsity_pattern::SparseMatrixCSC{Bool}
    jointtv_nzmatdiffs::Vector{SparseMatrixCSC}
end

get_propensity_gradients(SA::ForwardSensFspMatrixSparse) = SA.propensity_gradients
get_timeinvariant_matdiffs(SA::ForwardSensFspMatrixSparse) = SA.timeinvariant_matdiffs
get_separabletv_sparsity_pattern(SA::ForwardSensFspMatrixSparse) = SA.separabletv_matdiff_sparsity_pattern
get_separabletv_nzmatdiffs(SA::ForwardSensFspMatrixSparse) = SA.separabletv_nzmatdiffs
get_separabletv_mat_map(SA::ForwardSensFspMatrixSparse) = SA.separabletv_mat_map
get_jointtv_sparsity_pattern(SA::ForwardSensFspMatrixSparse) = SA.jointtv_matdiff_sparsity_pattern
get_jointtv_nzmatdiffs(SA::ForwardSensFspMatrixSparse) = SA.jointtv_nzmatdiffs

function ForwardSensFspMatrixSparse{RealT}(model::CmeModelWithSensitivity{IntT}, statespace::StateSpaceSparse{NS,NR,IntT,SizeT}) where {NS,NR,IntT<:Integer,SizeT<:Integer,RealT<:AbstractFloat}
    parameter_count = get_parameter_count(model)
    grad_sparsity_patterns = get_gradient_sparsity_patterns(model)
    propensity_gradients = get_propensity_gradients(model)
    parameters = get_parameters(model)

    fspmatrix = FspMatrixSparse{RealT}(statespace, get_propensities(model), parameters = get_parameters(model))
    timeinvariant_ids = get_timeinvariant_propensity_ids(fspmatrix)
    jointtv_ids = get_jointtv_propensity_ids(fspmatrix)
    separabletv_ids = get_separabletv_propensity_ids(fspmatrix)

    n = get_rowcount(fspmatrix)
    timeinvariant_matdiffs = []
    tinvariant_sparsity = grad_sparsity_patterns[timeinvariant_ids, :]
    rowvals = SparseArrays.getrowval(tinvariant_sparsity)
    for ip in 1:parameter_count
        rids = Vector{SizeT}()
        cids = Vector{SizeT}()
        vals = Vector{RealT}()
        for i in nzrange(tinvariant_sparsity, ip)
            r = timeinvariant_ids[rowvals[i]]
            Is, Js, Vs = _generate_sparsematrix_entries(statespace, get_single_pardiff(propensity_gradients[r], ip), r, parameters, valtype = RealT)
            append!(rids, Is)
            append!(cids, Js)
            append!(vals, Vs)
        end
        push!(timeinvariant_matdiffs, sparse(rids, cids, vals, n, n))
    end

    separabletv_sparsity = grad_sparsity_patterns[separabletv_ids, :]
    separabletv_nzmatdiffs = []
    separabletv_mat_map = []
    rowvals = SparseArrays.getrowval(separabletv_sparsity)
    for ip in 1:parameter_count
        for i in nzrange(separabletv_sparsity, ip)
            r = separabletv_ids[rowvals[i]]
            Is, Js, Vs = _generate_sparsematrix_entries(statespace, get_single_statefactor_pardiff(propensity_gradients[r], ip), r, parameters, valtype = RealT)
            push!(separabletv_nzmatdiffs, sparse(Is, Js, Vs, n, n))
            push!(separabletv_mat_map, findfirst(separabletv_ids .== r))
        end
    end

    jointtv_sparsity = grad_sparsity_patterns[jointtv_ids, :]
    jointtv_nzmatdiffs = []
    rowvals = SparseArrays.getrowval(jointtv_sparsity)
    for ip in 1:parameter_count
        for i in nzrange(jointtv_sparsity, ip)
            r = jointtv_ids[rowvals[i]]
            Is, Js, Vs = _generate_sparsematrix_entries(statespace, nothing, r, parameters, valtype = RealT)
            push!(jointtv_nzmatdiffs, sparse(Is, Js, Vs, n, n))
        end
    end


    return ForwardSensFspMatrixSparse{NS,NR,IntT,RealT}(
        fspmatrix,
        propensity_gradients,
        timeinvariant_matdiffs,
        separabletv_sparsity,
        separabletv_nzmatdiffs,
        separabletv_mat_map,
        jointtv_sparsity,
        jointtv_nzmatdiffs
    )
end

function matvec!(out::AbstractVector{RealT}, t::AbstractFloat, SA::ForwardSensFspMatrixSparse, vs::AbstractVector{RealT}) where {RealT<:AbstractFloat}
    A = SA.fspmatrix
    θ = get_parameters(A)
    parameter_count = length(θ)
    propensities = get_propensities(A)
    propensity_gradients = get_propensity_gradients(SA)

    jointtv_ids = get_jointtv_propensity_ids(A)
    separabletv_ids = get_separabletv_propensity_ids(A)

    ti_matdiffs = get_timeinvariant_matdiffs(SA)

    septv_pattern = get_separabletv_sparsity_pattern(SA)
    septv_map = get_separabletv_mat_map(SA)
    septv_matdiffs = get_separabletv_nzmatdiffs(SA)

    jointtv_pattern = get_jointtv_sparsity_pattern(SA)
    jointtv_matdiffs = get_jointtv_nzmatdiffs(SA)

    n = get_rowcount(SA.fspmatrix)
    matvec!(view(out, 1:n), t, A, view(vs, 1:n))
    for ip in 1:parameter_count
        outview = view(out, ip*n+1:(ip+1)*n)
        vsview = view(vs, ip*n+1:(ip+1)*n)
        matvecadd!(outview, t, A, vsview)

        mul!(outview, ti_matdiffs[ip], view(vs, 1:n))        

        for i in nzrange(septv_pattern, ip)
            r = separabletv_ids[SparseArrays.getrowval(septv_pattern)[i]]

            mul!(outview, septv_matdiffs[i], view(vs, 1:n), propensities[r].tfactor(t, θ), 1.0)
            
            j = septv_map[i]
            tfactordiff = get_single_tfactor_pardiff(propensity_gradients[r], ip)
            mul!(outview, A.separabletv_factormatrices[j], view(vs, 1:n), tfactordiff(t, θ), 1.0)
        end

        for i in nzrange(jointtv_pattern, ip)
            r = jointtv_ids[SparseArrays.getrowval(jointtv_pattern)[i]]
            pdiff = get_single_pardiff(propensity_gradients[r], ip)
            _update_sparsematrix!(jointtv_matdiffs[i], A.states, pdiff, t, θ)
            mul!(outview, jointtv_matdiffs[i], view(vs, 1:n), 1, 1)
        end
    end
    return nothing
end

