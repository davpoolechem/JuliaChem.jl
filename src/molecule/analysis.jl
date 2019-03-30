import Base.Threads

"""
     coordinate_analysis(coord::Array{Float64,2})
Summary
======
Perform the core molecular coordinate analysis algorithm.

Arguments
======
coord = molecular coordinates
"""
function coordinate_analysis(coord::Array{Float64,2})
    #determine some pre-information
    natoms::Int64 = size(coord)[1]

    #bond lengths
    bond_lengths::Array{Float64,2} = zeros(natoms,natoms)

    Threads.@threads for iatom in 1:natoms
        for jatom in 1:natoms
            diff::Array{Float64,1} = (coord[iatom,2:4].-coord[jatom,2:4]).^2
            bond_lengths[iatom,jatom] = sqrt(reduce(+,diff))
        end
    end

    println("----------------------------------------          ")
    println("        Printing bond lengths...                  ")
    println("----------------------------------------          ")
    println(" ")
    println("Atom #1   Atom #2     Bond length")
    for iatom in 1:natoms, jatom in 1:natoms
        if (iatom > jatom)
            println("   ",iatom,"         ",jatom,"     ",bond_lengths[iatom,jatom])
        end
    end
    println(" ")

    #bond angles

end
