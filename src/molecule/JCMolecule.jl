"""
     module JCMolecule
The module required for determination of molecular coordinate-based properties
(such as bond lengths, bond angles, and dihedral angles). Import this module
into the script when you wish to determine such information. Note that this
module is not strictly necessary for every calculation.
"""
module JCMolecule

Base.include(@__MODULE__,"MoleculeAnalysis.jl")

"""
     run(input_info::Dict{String,Dict{String,Any}})
Execute the JuliaChem molecular coordinate analysis functions.

One input variable is required:
1. input_info = Information gathered from the input file.

No variables are output.

Thus, proper use of the Molecule.run() function would look like this:

```
Molecule.run(input_info)
```
"""
function run(mol::MolStructs.Molecule)
    comm=MPI.COMM_WORLD

    if (MPI.Comm_rank(comm) == 0)
        println("---------------------------------------------------------------------------------")
        println("                       ========================================          ")
        println("                             MOLECULAR COORDINATE ANALYSIS               ")
        println("                       ========================================          ")
        println("")
    end

    #== print coordinates ==#
    print_xyz(mol) 

    #== compute bond lengths ==#
    bond_lengths = analyze_bond_lengths(mol)

    #== compute bond angles ==#
    #bond_angles = analyze_bond_angles(mol,bond_lengths)

    if (MPI.Comm_rank(comm) == 0)
        println("                       ========================================          ")
        println("                                END COORDINATE ANALYSIS                  ")
        println("                       ========================================          ")
    end
end
export run

end
