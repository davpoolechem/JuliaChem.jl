module Molecule

import MPI

using MoleculeAnalysis
"""
     do_coordinate_analysis(coord::Array{Float64,2})
Summary
======
Execute the JuliChem molecular coordinate analysis algorithm.

Arguments
======
coord = molecular coordinates
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
export do_coordinate_analysis

end
