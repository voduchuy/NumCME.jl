using Julifsp
using Test


@testset "CME model" begin
    include("test_propensity.jl")    
end

@testset "Sensitivity model" begin 
    include("test_autodiff.jl")
end

@testset "State Space" begin
    include("test_statespace.jl")    
end

@testset "Fsp Matrix" begin 
    include("test_fspmat.jl")
end

@testset "Transient CME" begin 
    include("test_solver.jl")
end

@testset "ForwradSensFspMatrix" begin 
    include("test_sensmat.jl")
end
