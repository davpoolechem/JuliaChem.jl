Base.include(@__MODULE__,"orbital_energies.jl")

"""
     properties(scf::Data,FLAGS::Flags)
Summary
======
Compute properties for RHF wave function.

Arguments
======
scf = Core HF data structures
"""
function do_properties(scf::Data,flags::Flags)
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
