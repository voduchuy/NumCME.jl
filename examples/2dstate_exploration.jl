using BenchmarkTools
using ChemicalMasterEquations 


function explore_states(state_type::Type{ST}, index_type::Type{T}) where {ST<:Integer, T<: Integer}
    S = Matrix{state_type}([[1;0] [-1;0] [0;1] [0;-1]])
    statespace = StateSpaceSparse(S, Vector{state_type}([0,0]); index_type=index_type) 
    expand!(statespace, 200)
    return statespace
end

@btime space16 = explore_states(Int32, UInt32);
@btime space64 = explore_states(Int64, UInt64);


