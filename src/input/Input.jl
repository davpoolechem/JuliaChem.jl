"""
     module Input
The module required for reading in and processing the selected input file.
Import this module into the script when you need to process an input file
(which will be every single calculation).
"""
module Input

using Base.Threads
using Distributed

using InputFunctions
using InputStructs

"""
     run()
Performs operations necessary to read in, process, and extract data from the
selected input file. The core function from the Input module; put Input.run()
into the script to read in the input file.

No input variables are required.

Two variables are output:
1. flags = The calculation flags from the input file.
2. coord = The molecular coordinates.

Thus, proper use of the Input.run() would look like this:
>flags, coord = Input.run()
"""
function run()
    #read in .inp and .dat files
    println("-------------------------------------------------------------------------------------")
    println("                       ========================================          ")
    println("                                READING INPUT DATA FILE                  ")
    println("                       ========================================          ")
    println(" ")

    directory::String = pwd()
    #println("Input file: ", directory*"/"*input_file)
    println(" ")
    println("Number of worker processes: ", Distributed.nworkers())
    println("Number of threads: ", Threads.nthreads())
    flags::Flags = input_flags()

    coord::Array{Float64,2} = input_coord()

    println(" ")
    println("                       ========================================          ")
    println("                                       END INPUT                         ")
    println("                       ========================================          ")
    return (flags,coord)
end
export do_input

end
