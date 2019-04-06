module Properties

using OrbitalEnergies

using InputStructs
using RHFStructs

import MPI

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
    comm=MPI.COMM_WORLD

    if (MPI.Comm_rank(comm) == 0)
        println("--------------------------------------------------------------------------------------")
        println("                       ========================================          ")
        println("                                   SYSTEM PROPERTIES                     ")
        println("                       ========================================          ")
        println("")
    end

    orbital_energies(scf,flags)

    if (MPI.Comm_rank(comm) == 0)
        println("                       ========================================          ")
        println("                                 END SYSTEM PROPERTIES                   ")
        println("                       ========================================          ")
    end
end
export do_properties

end
