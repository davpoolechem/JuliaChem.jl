"""
     module JCMolecule
The module required for determination of molecular coordinate-based properties
(such as bond lengths, bond angles, and dihedral angles). Import this module
into the script when you wish to determine such information. Note that this
module is not strictly necessary for every calculation.
"""
module JCMolecule

using MoleculeAnalysis

import MPI
"""
     run(coord::Array{Float64,2})
Execute the JuliaChem molecular coordinate analysis functions.

One input variable is required:
1. coord = The molecular coordinates.

No variables are output.

Thus, proper use of the Molecule.run() function would look like this:
>Molecule.run(coord)
"""
function run(coord::Array{Float64,2})
    comm=MPI.COMM_WORLD

    if (MPI.Comm_rank(comm) == 0)
        println("--------------------------------------------------------------------------------------")
        println("                       ========================================          ")
        println("                             MOLECULAR COORDINATE ANALYSIS               ")
        println("                       ========================================          ")
        println("")
    end

    coordinate_analysis(coord)

    if (MPI.Comm_rank(comm) == 0)
        println("                       ========================================          ")
        println("                                END COORDINATE ANALYSIS                  ")
        println("                       ========================================          ")
    end
end
export run

end
