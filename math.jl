import LinearAlgebra

function ∑(array::Array{Float64})
    return sum(array)
end

function ∑(array_1::Array{Float64},array_2::Array{Float64})
    return LinearAlgebra.dot(array_1,array_2)
end
