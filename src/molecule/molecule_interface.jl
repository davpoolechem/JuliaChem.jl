Base.include(@__MODULE__, "analysis.jl")

"""
     do_coordinate_analysis(coord::Array{Float64,2})
Summary
======
Execute the JuliChem molecular coordinate analysis algorithm.

Arguments
======
coord = molecular coordinates
"""
function do_coordinate_analysis(coord::Array{Float64,2})
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