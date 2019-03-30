#------------------------------------#
#            integrals.jl            #
#------------------------------------#
"""
    read_in_oei(data::Array{String,1}, type::String)
Summary
======
Extract one-electron integrals from data file object. Type dictates which
one-electron integrals are returned.

Arguments
======
data = name of data file object to process

type = tag for determining which type of one-electron integrals are returned:
1. type=OVR returns overlap integrals.
2. type=KEI returns kinetic energy integrals.
3. type=NAI returns nuclear-electron attraction integrals.
"""
function read_in_oei(oei::Array{Float64,2})
    nbf::Int64 = 7
    nbf2::Int64 = nbf*(nbf+1)/2

    oei_matrix::Array{Float64,2} = zeros(nbf,nbf)
    for index::Int64 in 1:(nbf2/2)
        i::Int64 = oei[index,1]
        j::Int64 = oei[index,2]

        oei_matrix[i,j] = oei[index,3]
        oei_matrix[j,i] = oei[i,j]
    end

    return oei_matrix
end

"""
    read_in_tei(data::Array{String,1})
Summary
======
Extract two-electron integrals from data file object.

Arguments
======
data = name of data file object to process
"""
function read_in_tei(tei::Array{Float64,2})
    nint::Int64 = 228

    tei_array::Array{Float64,1} = zeros(2401)
    for index::Int64 in 1:nint
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
