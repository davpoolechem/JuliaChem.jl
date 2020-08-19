"""
     module JCProperties
The module required for computation of a variety of properties, including
dipole moment, Mulliken charges, and orbital energies. Import this module
into the script when you wish to determine such information. Note that this
module is not strictly necessary for every calculation.
"""
module JCProperties

#Base.include(@__MODULE__,"OrbitalEnergies.jl")

using MPI

"""
     run(scf::Data,input_info::Dict{String,Dict{String,Any}})
Compute the dipole moment, Mulliken charges, and orbital energies of the
system in question.

Two input variables are required:
1. scf = Data saved from the SCF calculation.
2. input_info = Information gathered from the input file.

No variables are output.

Thus, proper use of the Properties.run() function would look like this:

```
Properties.run(scf, input_info)
```
"""
function run(scf::Data,input_info::Dict{String,Dict{String,Any}})
    comm=MPI.COMM_WORLD

    if (MPI.Comm_rank(comm) == 0)
        println("--------------------------------------------------------------------------------------")
        println("                       ========================================          ")
        println("                                   SYSTEM PROPERTIES                     ")
        println("                       ========================================          ")
        println("")
    end

    #set up rhf flags
    basis_info::Dict{String,Any} = input_info["Basis Flags"]
    basis_flags::Basis_Flags = Basis_Flags(basis_info["norb"],
      basis_info["nocc"])

    orbital_energies(scf,basis_flags)

    if (MPI.Comm_rank(comm) == 0)
        println("                       ========================================          ")
        println("                                 END SYSTEM PROPERTIES                   ")
        println("                       ========================================          ")
    end
end
export run

end
