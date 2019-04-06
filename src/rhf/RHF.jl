module RHF

import MPI

using RHFSCF
using RHFStructs

using InputStructs

"""
    do_rhf(dat::Array{String,1})
Summary
======
Execute the JuliChem RHF algorithm.

Arguments
======
dat = input data file object
"""
function run(flags::Flags)
    comm=MPI.COMM_WORLD

    if (MPI.Comm_rank(comm) == 0)
        println("--------------------------------------------------------------------------------------")
        println("                       ========================================          ")
        println("                         RESTRICTED CLOSED-SHELL HARTREE-FOCK            ")
        println("                       ========================================          ")
        println("")
    end

    #GC.enable(false)
    scf::Data = rhf_energy(flags)
    #GC.enable(true)
    #GC.gc()

    if (MPI.Comm_rank(comm) == 0)
        println("                       ========================================          ")
        println("                             END RESTRICTED CLOSED-SHELL                 ")
        println("                                     HARTREE-FOCK                        ")
        println("                       ========================================          ")
    end

    return scf
end
export do_rhf

end
