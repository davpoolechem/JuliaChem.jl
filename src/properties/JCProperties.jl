"""
     module JCProperties
The module required for computation of a variety of properties, including
dipole moment, Mulliken charges, and orbital energies. Import this module
into the script when you wish to determine such information. Note that this
module is not strictly necessary for every calculation.
"""
module JCProperties

using OrbitalEnergies

using InputStructs
using RHFStructs

import MPI

"""
     run(scf::Data,flags::Flags)
Compute the dipole moment, Mulliken charges, and orbital energies of the
system in question.

Two input variables are required:
1. scf = Data saved from the SCF calculation.
2. flags = The calculation flags from the input file.

No variables are output.

Thus, proper use of the Properties.run() function would look like this:

```
Properties.run(scf, flags)
```
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
export run

end
