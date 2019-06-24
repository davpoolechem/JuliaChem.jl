#------------------------------#
#            Math.jl           #
#------------------------------#
using LinearAlgebra

#=
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
=#
function ∑(array::Array{T}) where {T<:Number}
  return sum(array)
end
export ∑

#=
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
=#
function ∑(array_1::Array{T,1},array_2::Array{T,1}) where {T<:Number}
  return LinearAlgebra.dot(array_1,array_2)
end
#=
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
=#
function ∑(array_1::Array{T,2},array_2::Array{T,2}) where {T<:Number}
  return LinearAlgebra.dot(array_1,array_2)
end
export ∑

#=
"""
    Method #4 of ∑:
    ∑(array_1::Array{Float64,2},array_2::Array{Float64,2},array_center::Array{Float64,2})
Summary
======
Return quadratic transform of three two-dimensional arrays.

Arguments
======
array_1 = first array in the dot product

array_2 = second array in the dot product

array_center = center array in the dot product

"""
=#
function ∑(array_1::Array{T},array_2::Array{T},array_center::Array{T,2}) where {T<:Number}
  array::Float64 = transpose(array_2)*array_center*array_1
  #array = (array .* array).^0.5
  #return reduce(+,array)
  return array
end
