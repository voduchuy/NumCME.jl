export AbstractFspMatrix, matvec, matvec!

"""
Abstract type for the FSP-truncated infinitesimal generator of the Chemical Master Equation.
"""
abstract type AbstractFspMatrix end



include("sparse/fspsparsematrix.jl")


