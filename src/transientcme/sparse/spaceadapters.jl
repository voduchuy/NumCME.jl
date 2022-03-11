using DifferentialEquations.DiffEqBase: DEIntegrator

export AbstractSpaceAdapterSparse, AbstractSpaceAdapter, init!, adapt!

abstract type AbstractSpaceAdapter end
abstract type AbstractSpaceAdapterSparse <: AbstractSpaceAdapter end 

include("rstepadapters.jl")


