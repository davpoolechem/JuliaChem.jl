Base.include(@__MODULE__,"../math/math.jl")

using JCStructs

using MPI

function orbital_energies(scf::Data,FLAGS::Basis_Flags)
    norb = FLAGS.NORB
    comm = MPI.COMM_WORLD

    F_mo::Array{Float64,2} = zeros(norb,norb)
    for i::Int64 in 1:norb, j::Int64 in 1:i
        #F_mo[i,j] = ∑(scf.Coeff[:,1:j],scf.Coeff[:,1:i],scf.Fock)
        F_mo[i,j] = ∑(scf.Coeff[:,j],scf.Coeff[:,i],scf.Fock)
        F_mo[j,i] = F_mo[i,j]
    end
    if (MPI.Comm_rank(comm) == 0)
        println("----------------------------------------          ")
        println(" Printing molecular orbital energies...           ")
        println("----------------------------------------          ")
        println(" ")
        println("Orbital #     Orbital energy")
        for index::Int64 in 1:norb
            println("   ",index,"       ",F_mo[index,index])
        end
        println(" ")
    end
end
