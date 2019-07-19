Base.include(@__MODULE__,"../basis/BasisStructs.jl")
using .BasisStructs

using Base.Threads

function read_in_enuc()
	enuc::Float64 = input_enuc()

	return enuc
end
#=
"""
		get_oei_matrix(oei::Array{Float64,2})
Summary
======
Extract one-electron integrals from data file object. Kinetic energy integrals,
overlap integrals, and nuclear attraction integrals can all be extracted.

Arguments
======
oei = array of one-electron integrals to extract
"""
=#
function read_in_oei(oei::Array{Any,1}, nbf::Int64)
	nbf2::Int64 = nbf*(nbf+1)/2

    ioff::Array{Int64,1} = map((x) -> x*(x-1)/2, collect(1:nbf*(nbf+1)))

	oei_matrix::Array{Float64,2} = Matrix{Float64}(undef,(nbf,nbf))
	Threads.@threads for index::Int64 in 1:nbf2
        i::Int64 = ceil(((-1+sqrt(1+8*index))/2))
        j::Int64 = index - ioff[i]

        eri = oei[index]
		oei_matrix[i,j] = oei[index]
		oei_matrix[j,i] = oei_matrix[i,j]
	end

	return oei_matrix
end
