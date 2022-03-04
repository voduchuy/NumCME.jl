using DifferentialEquations.DiffEqBase: DEIntegrator

export SparseSpaceAdapter, SpaceAdapter, init!, adapt!

abstract type SpaceAdapter end
abstract type SparseSpaceAdapter <: SpaceAdapter end 

include("rstepadapters.jl")


