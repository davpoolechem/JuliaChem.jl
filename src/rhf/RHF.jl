"""
     module Properties
The module required for computation of the wave function using the *Restricted
Hartree-Fock* (RHF) method in a Self-Consistent Field (SCF) calculation. This
module will be used often, as the RHF wave function is often the zeroth-order
wave function for closed-shell systems.
"""
module RHF

using RHFSCF
using RHFStructs

using InputStructs

"""
    run(flags::Flags)

Execute the JuliChem RHF algorithm.

One input variable is required:
1. flags = The calculation flags from the input file.

One variable is output:
1. scf = Data saved from the SCF calculation.

Thus, proper use of the RHF.run() function would look like this:
>scf = RHF.run(flags)
"""
function run(flags::Flags)
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
