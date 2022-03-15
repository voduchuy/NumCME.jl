using NumCME 
using StaticArrays 
using Test 

states = [[@MVector([0, i]) for i in 0:10];[@MVector([1, i]) for i in 0:10]]
vals = [[0.0 for i in 0:10]; [1.0 for i in 0:10]]
v = FspVectorSparse(states, vals)
v1 = sum(v, [1])
v2 = sum(v, [2])
@test typeof(v) <: AbstractFspVector
@test typeof(v) <: FspVectorSparse{2, <:Integer, <:Real}
@test sum(v) == 11.0
@test nnz(v1) == 11
@test nnz(v2) == 2

@test_throws ArgumentError FspVectorSparse([@MVector([0, 1])], [1.0,2.0])
@test_throws ArgumentError FspVectorSparse([@MVector([0,1]), @MVector([1,2])], [1.0])