using EtFsp
using Test

@testset "StateSpace" begin
    include("test_statespace.jl")    
end

@testset "FspMatrix" begin 
    include("test_fspmat.jl")
end
