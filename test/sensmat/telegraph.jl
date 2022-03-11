using Julifsp
using Test
using BenchmarkTools
using SparseArrays
using LinearAlgebra: mul!, norm

𝕊 = [[-1, 1, 0] [1, -1, 0] [0, 0, 1] [0, 0, -1]]

x₀ = [1, 0, 0]
k₀₁ = 0.05
k₁₀ = 0.1
λ = 5.0
γ = 0.5
L = 20.0
θ = [k₀₁, k₁₀, λ, γ, L]
a1 = propensity() do x, p
    p[1] * x[1]
end
a2 = propensity() do x, p
    p[2] * x[2]
end
a3 = propensity() do x, p
    p[3] * x[2]
end
a4 = propensity((t, p) -> max(0.0, 1.0 - sin(π * t / p[5]))) do x, p
    p[4] * x[3]
end
propensities = [a1, a2, a3, a4]
RNACOUNT_MAX = 500
sensmodel = CmeModelWithSensitivity(CmeModel(𝕊, propensities, θ))
test_space = StateSpaceSparse(get_stoich_matrix(sensmodel), [[[1, 0, i] for i in 0:RNACOUNT_MAX]; [[0, 1, i] for i in 0:RNACOUNT_MAX]])

## Test that output vectors of mavec sum to ≈ 0
function test_zero_sum()
    pcount = get_parameter_count(sensmodel)
    sensmat = ForwardSensFspMatrixSparse{Float64}(sensmodel, test_space);
    n = get_rowcount(sensmat.fspmatrix)
    v = ones(n + pcount * n)
    out = similar(v)
    matvec!(out, 20.0, sensmat, v)
    return isapprox(sum(out) / sum(v), 0.0, atol = n*(pcount+1)*eps())
end
@test test_zero_sum()

## Test against analytic matrices
function getA(nmax::Integer, t::AbstractFloat, parameters::AbstractVector{<:AbstractFloat})    
    rowids = Vector{UInt32}()
    colids = Vector{UInt32}()
    nzvals = Vector{Float64}()
    for gene in 0:1
        for rnacount in 0:nmax 
            x = [1-gene, gene, rnacount]
            idx = gene*(nmax+1) + rnacount+1 
            # Diagonal
            push!(rowids, idx)
            push!(colids, idx)
            push!(nzvals, -1.0*(a1(x, parameters) + a2(x, parameters) + a3(x, parameters) + a4(t, x, parameters)))
            # Offdiagonal for goff->gon 
            if gene == 0
                push!(colids, idx)
                push!(rowids, nmax+1 + rnacount+1)
                push!(nzvals, a1(x, parameters))
            end
            # Offdiagonal for gon->goff 
            if gene == 1
                push!(colids, idx)
                push!(rowids, rnacount+1)
                push!(nzvals, a2(x, parameters))
            end
            # Offdiagonal for gon->gon+rna 
            if gene == 1
                push!(colids, idx)
                push!(rowids, (rnacount < nmax) ? nmax+1+rnacount+2 : 2*(nmax+1)+3)
                push!(nzvals, a3(x, parameters))
            end
            # Offdiagonal for rna->0
            if rnacount > 0
                push!(colids, idx)
                push!(rowids, gene*(nmax+1) + rnacount)
                push!(nzvals, a4(t, x, parameters))
            end
        end
    end
    return sparse(rowids, colids, nzvals, 2*(nmax+1)+4, 2*(nmax+1)+4)
end
# Generate the partial derivative matrix wrt k10 corresponding to a rectangular state space 
function ∂A∂k01(nmax::Integer)
    sink_count = 4
    state_count = 2 * (nmax + 1)
    rowids = Vector{UInt32}()
    colids = Vector{UInt32}()
    nzvals = Vector{Float64}()
    for rnacount in 0:nmax
        # Fill column corresponding to state (1, 0, rnacount)
        push!(rowids, rnacount + 1)
        push!(colids, rnacount + 1)
        push!(nzvals, -1.0)
        push!(rowids, nmax + 1 + rnacount + 1)
        push!(colids, rnacount + 1)
        push!(nzvals, 1.0)
    end
    return sparse(rowids, colids, nzvals, state_count + sink_count, state_count + sink_count)
end
function ∂A∂k10(nmax::Integer)
    sink_count = 4
    state_count = 2 * (nmax + 1)
    rowids = Vector{UInt32}()
    colids = Vector{UInt32}()
    nzvals = Vector{Float64}()
    for rnacount in 0:nmax
        # Fill column corresponding to state (0, 1, rnacount)
        push!(colids, nmax + 1 + rnacount + 1)
        push!(rowids, nmax + 1 + rnacount + 1)
        push!(nzvals, -1.0)
        push!(colids, nmax + 1 + rnacount + 1)
        push!(rowids, rnacount + 1)
        push!(nzvals, 1.0)
    end
    return sparse(rowids, colids, nzvals, state_count + sink_count, state_count + sink_count)
end
function ∂A∂λ(nmax::Integer)
    sink_count = 4
    state_count = 2 * (nmax + 1)
    rowids = Vector{UInt32}()
    colids = Vector{UInt32}()
    nzvals = Vector{Float64}()
    for rnacount in 0:nmax
        # Fill column corresponding to state (0, 1, rnacount)
        push!(colids, nmax + 1 + rnacount + 1)
        push!(rowids, nmax + 1 + rnacount + 1)
        push!(nzvals, -1.0)
        push!(colids, nmax + 1 + rnacount + 1)
        push!(rowids, (rnacount < nmax) ? nmax + 1 + rnacount + 2 : state_count + 3)
        push!(nzvals, 1.0)
    end
    return sparse(rowids, colids, nzvals, state_count + sink_count, state_count + sink_count)
end
function ∂A∂γ(nmax::Integer, t, p)
    sink_count = 4
    state_count = 2 * (nmax + 1)
    rowids = Vector{UInt32}()
    colids = Vector{UInt32}()
    nzvals = Vector{Float64}()
    for rnacount in 1:nmax
        # Fill column corresponding to state (1, 0, rnacount)
        push!(rowids, rnacount + 1)
        push!(colids, rnacount + 1)
        push!(nzvals, -max(0.0, 1.0 - sin(π * t / p[5])) * rnacount)
        push!(rowids, rnacount)
        push!(colids, rnacount + 1)
        push!(nzvals, max(0.0, 1.0 - sin(π * t / p[5])) * rnacount)

        # Fill column corresponding to state (0, 1, rnacount)
        push!(colids, nmax + 1 + rnacount + 1)
        push!(rowids, nmax + 1 + rnacount + 1)
        push!(nzvals, -max(0.0, 1.0 - sin(π * t / p[5])) * rnacount)
        push!(rowids, nmax + 1 + rnacount)
        push!(colids, nmax + 1 + rnacount + 1)
        push!(nzvals, max(0.0, 1.0 - sin(π * t / p[5])) * rnacount)

    end
    return sparse(rowids, colids, nzvals, state_count + sink_count, state_count + sink_count)
end
function ∂A∂L(nmax::Integer, t, p)
    sink_count = 4
    state_count = 2 * (nmax + 1)
    rowids = Vector{UInt32}()
    colids = Vector{UInt32}()
    nzvals = Vector{Float64}()
    for rnacount in 1:nmax
        # Fill column corresponding to state (1, 0, rnacount)
        push!(rowids, rnacount + 1)
        push!(colids, rnacount + 1)
        push!(nzvals, (1.0-sin(π*t/p[5]) >= 0) ? -1.0*p[4]*rnacount*cos(π*t/p[5])*π*t/(p[5]^2) : 0.0)
        push!(rowids, rnacount)
        push!(colids, rnacount + 1)
        push!(nzvals, (1.0-sin(π*t/p[5]) >= 0) ? 1.0*p[4]*rnacount*cos(π*t/p[5])*π*t/(p[5]^2) : 0.0)

        # Fill column corresponding to state (0, 1, rnacount)
        push!(colids, nmax + 1 + rnacount + 1)
        push!(rowids, nmax + 1 + rnacount + 1)
        push!(nzvals, (1.0-sin(π*t/p[5]) >= 0) ? -1.0*p[4]*rnacount*cos(π*t/p[5])*π*t/(p[5]^2) : 0.0)
        push!(rowids, nmax + 1 + rnacount)
        push!(colids, nmax + 1 + rnacount + 1)
        push!(nzvals, (1.0-sin(π*t/p[5]) >= 0) ? 1.0*p[4]*rnacount*cos(π*t/p[5])*π*t/(p[5]^2) : 0.0)
    end
    return sparse(rowids, colids, nzvals, state_count + sink_count, state_count + sink_count)
end

function test_against_analytic_matrix(t_test::AbstractFloat, rnacount_max::Integer)
    sensmat = ForwardSensFspMatrixSparse{Float64}(sensmodel, test_space)
    pcount = get_parameter_count(sensmodel)
    n = get_rowcount(sensmat.fspmatrix)
    v = ones(n + pcount * n)
    v ./= sum(v)
    out_sensmat = similar(v)
    matvec!(out_sensmat, t_test, sensmat, v)

    A = getA(rnacount_max, t_test, θ)
    dAs = []
    push!(dAs, ∂A∂k01(rnacount_max))
    push!(dAs, ∂A∂k10(rnacount_max))
    push!(dAs, ∂A∂λ(rnacount_max))
    push!(dAs, ∂A∂γ(rnacount_max, t_test, θ))
    push!(dAs, ∂A∂L(rnacount_max, t_test, θ))
    out_ref = similar(v)
    mul!(view(out_ref, 1:n), A, view(v, 1:n))
    for ip in 1:pcount 
        mul!(view(out_ref, ip*n+1:(ip+1)*n), dAs[ip], view(v, 1:n))
        mul!(view(out_ref, ip*n+1:(ip+1)*n), A, view(v, ip*n+1:(ip+1)*n), 1, 1)
    end
    return isapprox(norm(out_sensmat- out_ref, Inf), 0.0, atol=n*eps())
end

for t_test in [10.0, 20.0, 30.0, 100.0]
    @test test_against_analytic_matrix(t_test, RNACOUNT_MAX)
end





