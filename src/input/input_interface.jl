Base.include(@__MODULE__, "../../input.jl")
import Base.Threads

#------------------------------#
#          JuliChem.jl         #
#------------------------------#
"""
     do_input(file::String)
Summary
======
Read in input file.

Arguments
======
file = name of input file to read in
"""
function do_input()
    #read in .inp and .dat files
    println("--------------------------------------------------------------------------------------")
    println("                       ========================================          ")
    println("                                READING INPUT DATA FILE                  ")
    println("                       ========================================          ")
    println(" ")

    directory::String = pwd()
    println("Input file: ", directory*"/"*input_file)
    println(" ")
    println("Number of worker processes: ", Distributed.nworkers())
    println("Number of threads: ", Threads.nthreads())
    flags::Flags = input_flags()

    coord::Array{Float64,2} = input_geometry()

    println(" ")
    println("                       ========================================          ")
    println("                                       END INPUT                         ")
    println("                       ========================================          ")
    return (flags,coord)
end
