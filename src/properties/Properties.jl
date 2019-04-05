"""
     module Properties
The module required for computation of a variety of properties, including
dipole moment, Mulliken charges, and orbital energies. Import this module
into the script when you wish to determine such information. Note that this
module is not strictly necessary for every calculation.
"""
module Properties

using OrbitalEnergies

using InputStructs
using RHFStructs

"""
     run(scf::Data,flags::Flags)
Compute the dipole moment, Mulliken charges, and orbital energies of the
system in question.

Two input variables are required:
1. scf = Data saved from the SCF calculation.
2. flags = The calculation flags from the input file.

No variables are output.

Thus, proper use of the Properties.run() function would look like this:
>Properties.run(scf, flags)
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
