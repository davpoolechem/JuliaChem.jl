import MPI
import Base.Threads
import LinearAlgebra

"""
     coordinate_analysis(coord::Array{Float64,2})
Summary
======
Perform the core molecular coordinate analysis algorithm.

Arguments
======
coord = molecular coordinates
"""
function analyze_bond_lengths(coord::Array{Float64,2})
    #determine some pre-information
    natoms::Int64 = size(coord)[1]
    comm=MPI.COMM_WORLD

    #calculate bond lengths
    bond_lengths::Array{Float64,2} = zeros(natoms,natoms)

    Threads.@threads for ijatom in 1:natoms*natoms
        iatom::Int64 = ceil(ijatom/natoms)
        jatom::Int64 =ijatom%natoms + 1

        if (iatom > jatom)
            diff::Array{Float64,1} = (coord[iatom,2:4].-coord[jatom,2:4]).^2
            bond_lengths[iatom,jatom] = sqrt(reduce(+,diff))
        end
    end

    #print bond lengths
    if (MPI.Comm_rank(comm) == 0)
        println("----------------------------------------          ")
        println("        Printing bond lengths...                  ")
        println("----------------------------------------          ")
        println(" ")
        println("Atom #1   Atom #2     Bond length")
    end
    for iatom in 1:natoms, jatom in 1:natoms
        if (iatom > jatom)
            if (MPI.Comm_rank(comm) == 0)
                println("   ",iatom,"         ",jatom,"     ",bond_lengths[iatom,jatom])
            end
        end
    end
    if (MPI.Comm_rank(comm) == 0)
        println(" ")
    end

    return bond_lengths
end

"""
     coordinate_analysis(coord::Array{Float64,2})
Summary
======
Perform the core molecular coordinate analysis algorithm.

Arguments
======
coord = molecular coordinates
"""
function analyze_bond_angles(coord::Array{Float64,2}, bond_lengths::Array{Float64,2})
    #determine some pre-information
    natoms::Int64 = size(coord)[1]
    comm=MPI.COMM_WORLD

    #calculate bond angles
    bond_angles::Array{Float64,3} = zeros(natoms,natoms,natoms)

    Threads.@threads for ijkatom in 1:natoms*natoms*natoms
        iatom::Int64 = ceil(ijkatom/(natoms*natoms))
        jatom::Int64 = ceil(ijkatom/natoms)%natoms+1
        katom::Int64 = (ijkatom%(natoms*natoms))%natoms + 1

        if (iatom > jatom && jatom > katom)
            if (bond_lengths[iatom,jatom] < 4.0 && bond_lengths[jatom,katom] < 4.0)
                e_ji = [ -(coord[jatom,1]-coord[iatom,1])/bond_lengths[iatom,jatom];
                         -(coord[jatom,2]-coord[iatom,2])/bond_lengths[iatom,jatom];
                         -(coord[jatom,3]-coord[iatom,3])/bond_lengths[iatom,jatom]; ]

                e_jk = [ -(coord[jatom,1]-coord[katom,1])/bond_lengths[jatom,katom];
                         -(coord[jatom,2]-coord[katom,2])/bond_lengths[jatom,katom];
                         -(coord[jatom,3]-coord[katom,3])/bond_lengths[jatom,katom]; ]

                e_ji = e_ji/LinearAlgebra.norm(e_ji)
                e_jk = e_jk/LinearAlgebra.norm(e_jk)

                bond_angles[iatom,jatom,katom] = (360/(2π))*acos(LinearAlgebra.dot(e_ji,e_jk))
            end
        end
    end

    #print bond angles
    if (MPI.Comm_rank(comm) == 0)
        println("----------------------------------------          ")
        println("         Printing bond angles...                  ")
        println("----------------------------------------          ")
        println(" ")
        println("Atom #1   Atom #2   Atom #3     Bond angle")
    end
    for iatom in 1:natoms, jatom in 1:natoms, katom in 1:natoms
        if (iatom > jatom && jatom > katom)
            if (bond_lengths[iatom,jatom] < 4.0 && bond_lengths[jatom,katom] < 4.0)
                if (MPI.Comm_rank(comm) == 0)
                    println("   ",iatom,"         ",jatom,"         ",katom,"     ",bond_angles[iatom,jatom,katom])
                end
            end
        end
    end
    if (MPI.Comm_rank(comm) == 0)
        println(" ")
    end

    return bond_lengths
end

function test(ijkatom::Int64)
    natoms = 3

    iatom::Int64 = ceil(ijkatom/(natoms*natoms))
    jatom::Int64 = ceil(ijkatom/natoms)%natoms+1
    katom::Int64 = (ijkatom%(natoms*natoms))%natoms + 1
    #println(iatom,", ",jatom,", ",katom)
    return (iatom, jatom, katom)
end

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
    #bond lengths
    bond_lengths::Array{Float64,2} = analyze_bond_lengths(coord)

    #bond angles
    analyze_bond_angles(coord,bond_lengths)
end
