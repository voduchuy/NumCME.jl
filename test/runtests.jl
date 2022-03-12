using Julifsp
using Test


@testset "CME model" begin
    include("test_propensity.jl")    
    include("test_autodiff.jl")
end

@testset "State Space" begin
    include("test_statespace.jl")    
end

@testset "Transient CME" begin 
    include("test_fspmat.jl")
    include("test_solver.jl")
end

@testset "Forward Sensitivity" begin 
    include("test_sensmat.jl")
    include("test_sensfsp.jl")
end
