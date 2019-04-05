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
system in question. The core function from the Properties module;
put Properties.run() into the script to compute aforementioned molecular
properties.
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
