module InputIntegrals

using Base.Threads

#------------------------------------#
#            integrals.jl            #
#------------------------------------#
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
function get_oei_matrix(oei::Array{Float64,2})
    nbf::Int64 = 1
    nbf2::Int64 = nbf*(nbf+1)/2

    oei_matrix::Array{Float64,2} = zeros(nbf,nbf)
    Threads.@threads for index::Int64 in 1:nbf2
        i::Int64 = oei[index,1]
        j::Int64 = oei[index,2]

        oei_matrix[i,j] = oei[index,3]
        oei_matrix[j,i] = oei_matrix[i,j]
    end

    return oei_matrix
end
export get_oei_matrix

#=
"""
    read_in_tei(data::Array{String,1})
Summary
======
Extract two-electron integrals from data file object.

Arguments
======
data = name of data file object to process
"""
=#
function get_tei_matrix(tei::Array{Float64,2})
    nint::Int64 = 1

    tei_array::Array{Float64,1} = zeros(2401)
    Threads.@threads for index::Int64 in 1:nint
        i::Int64 = tei[index,1]
        j::Int64 = tei[index,2]
        k::Int64 = tei[index,3]
        l::Int64 = tei[index,4]

        ij::Int64 = (i > j) ? i*(i+1)/2 + j : j*(j+1)/2 + i
        kl::Int64 = (k > l) ? k*(k+1)/2 + l : l*(l+1)/2 + k
        ijkl::Int64 = (ij > kl) ? ij*(ij+1)/2 + kl : kl*(kl+1)/2 + ij

        tei_array[ijkl] = tei[index,5]
    end

    return tei_array
end
export get_tei_matrix

end