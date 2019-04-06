"""
     module Molecule
The module required for determination of molecular coordinate-based properties
(such as bond lengths, bond angles, and dihedral angles). Import this module
into the script when you wish to determine such information. Note that this
module is not strictly necessary for every calculation.
"""
module Molecule

using MoleculeAnalysis

"""
     run(coord::Array{Float64,2})
Execute the JuliChem molecular coordinate analysis functions.

One input variable is required:
1. coord = The molecular coordinates.

No variables are output.

Thus, proper use of the Molecule.run() function would look like this:
>Molecule.run(coord)
"""
function run(coord::Array{Float64,2})
    println("--------------------------------------------------------------------------------------")
    println("                       ========================================          ")
    println("                             MOLECULAR COORDINATE ANALYSIS               ")
    println("                       ========================================          ")
    println("")

    coordinate_analysis(coord)

    println("                       ========================================          ")
    println("                                END COORDINATE ANALYSIS                  ")
    println("                       ========================================          ")
end
export do_coordinate_analysis

end
