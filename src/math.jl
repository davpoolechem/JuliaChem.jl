#------------------------------#
#            Math.jl           #
#------------------------------#
import LinearAlgebra

"""
    Method #1 of ∑:
    ∑(array::Array{Float64})
Summary
======
Return sum of elements in an array.

Arguments
======
array = array to be summed across
"""
function ∑(array::Array{Float64})
    return sum(array)
end

"""
    Method #2 of ∑:
    ∑(array_1::Array{Float64,1},array_2::Array{Float64,1})
Summary
======
Return dot product of two one-dimensional arrays.

Arguments
======
array_1 = first array in the dot product

array_2 = second array in the dot product

"""
function ∑(array_1::Array{Float64,1},array_2::Array{Float64,1})
    return LinearAlgebra.dot(array_1,array_2)
end

"""
    Method #3 of ∑:
    ∑(array_1::Array{Float64,2},array_2::Array{Float64,2})
Summary
======
Return dot product of two two-dimensional arrays.

Arguments
======
array_1 = first array in the dot product

array_2 = second array in the dot product

"""
function ∑(array_1::Array{Float64,2},array_2::Array{Float64,2})
    return LinearAlgebra.dot(array_1,array_2)
end
