module RHFInterface

using RHFSCF

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
function do_rhf(flags::Flags)
    println("--------------------------------------------------------------------------------------")
    println("                       ========================================          ")
    println("                         RESTRICTED CLOSED-SHELL HARTREE-FOCK            ")
    println("                       ========================================          ")
    println("")

    #GC.enable(false)
    scf::Data = rhf_energy(flags)
    #GC.enable(true)
    #GC.gc()

    println("                       ========================================          ")
    println("                             END RESTRICTED CLOSED-SHELL                 ")
    println("                                     HARTREE-FOCK                        ")
    println("                       ========================================          ")

    return scf
end
export do_rhf

end
