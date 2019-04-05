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
     do_coordinate_analysis(coord::Array{Float64,2})
Execute the JuliChem molecular coordinate analysis functions.
The core function from the Molecule module; put Molecule.run()
into the script to perform a coordinate analysis of the molecule.
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
