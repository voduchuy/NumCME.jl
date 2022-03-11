export AbstractStateSpace

"""
Abstract type for FSP state space. This is the supertype of all concrete FSP state space implementations.
"""
abstract type AbstractStateSpace end

include("sparse/sparsestatespace.jl")
