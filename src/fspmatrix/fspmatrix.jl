module fspmatrix     
    using EtFsp.statespace
    using SparseArrays
    import Base: *, size

    export Propensity, FspMatrix

    mutable struct Propensity         
        tfactor::Union{Function,Nothing}
        statefactor::Function         
    end

    mutable struct FspMatrix
        space::FspStateSpace
        rowcount::Int 
        colcount::Int

        propensity_functions::Vector{Propensity}
        tvarying_reactions::Vector{Int}
        
        tinvariant_matrix::Any
        tvarying_matrices::Vector{Any}
    end

    function FspMatrix(space::FspStateSpace{NS, NR, T}, propensity_functions::Vector{Propensity}) where {NS, NR, T<:Int}
        state_count = get_state_count(space)
        sink_count = get_sink_count(space)

        tvarying_reactions = Vector{Int}()
        tinvariant_reactions = Vector{Int}()
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

        rowindices = Vector{Int}()
        colindices = Vector{Int}()
        vals = Vector{Float64}()
        for reactionidx in tinvariant_reactions
            ridxs, cidxs, vs = _generate_sparsematrix_entries(space, propensity_functions[reactionidx].statefactor, reactionidx)
            rowindices = vcat(rowindices, ridxs)
            colindices = vcat(colindices, cidxs)
            vals = vcat(vals, vs)            
        end
        tinvariant_matrix = !isempty(tinvariant_reactions) ? sparse(rowindices, colindices, vals, n, n) : nothing 

        return FspMatrix(
            space,
            n,
            n,
            propensity_functions,
            tvarying_reactions,
            tinvariant_matrix,
            tvarying_matrices
        )
    end

    function _generate_sparsematrix_entries(space::FspStateSpace, statefactor::Function, reactionidx::Int)
        state_count = get_state_count(space)
        sink_count = get_sink_count(space)

        rowindices = Vector{Integer}()
        colindices = Vector{Integer}()
        vals = Vector{Float64}()

        sizehint!(rowindices, 2*state_count)
        sizehint!(colindices, 2*state_count)
        sizehint!(vals, 2*state_count)
        
        for idx=1:state_count

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
            push!(vals, -1.0*dval)
            
            # Populate the row corresponding to a reachable sink state (if any)
            sink_idx = space.sink_connectivity[idx][reactionidx]
            if sink_idx ≠ 0
                push!(rowindices, state_count + sink_idx)
                push!(colindices, idx)
                push!(vals, dval)
            end
        end                
        return rowindices[colindices .≠ 0], colindices[colindices .≠ 0], vals[colindices .≠ 0]
    end

    # Functions to query the dimensions of the matrix 
    function size(A::FspMatrix)
        return (A.rowcount, A.colcount)
    end

    function size(A::FspMatrix, dim::Int)
        if !(1≤ dim ≤ 2)
            throw(ArgumentError("Second argument must be either 1 or 2."))
        end
        return (dim == 1) ? A.rowcount : A.colcount
    end        

    # Matrix-vector multiplications
    function matmult_(t::AbstractFloat, A::FspMatrix, v::Vector{T}) where {T}        
        w = !(A.tinvariant_matrix isa Nothing) ? A.tinvariant_matrix*v : zeros(T, length(v))
        for (i, r) in enumerate(A.tvarying_reactions)
            w += A.propensity_functions[r].tfactor(t)*A.tvarying_matrices[i]*v
        end
        return w
    end

    function (A::FspMatrix)(t::AbstractFloat)
        return (t, A)
    end

    function *(A_at_t::Tuple{AbstractFloat, FspMatrix}, v::Vector)
        t = A_at_t[1]
        A = A_at_t[2]
        return matmult_(t, A, v)
    end

    function *(A::FspMatrix, v::Vector)
        return matmult_(0.0, A, v)
    end
end