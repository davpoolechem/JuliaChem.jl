module Properties

using OrbitalEnergies

using InputStructs
using RHFStructs

"""
     properties(scf::Data,FLAGS::Flags)
Summary
======
Compute properties for RHF wave function.

Arguments
======
scf = Core HF data structures
"""
function run(scf::Data,flags::Flags)
    println("--------------------------------------------------------------------------------------")
    println("                       ========================================          ")
    println("                                   SYSTEM PROPERTIES                     ")
    println("                       ========================================          ")
    println("")

    orbital_energies(scf,flags)

    println("                       ========================================          ")
    println("                                 END SYSTEM PROPERTIES                   ")
    println("                       ========================================          ")
end
export do_properties

end
