"""
     module JCInput
The module required for reading in and processing the selected input file.
Import this module into the script when you need to process an input file
(which will be every single calculation).
"""
module JCInput

include("InputFunctions.jl")
include("InputStructs.jl")

import MPI
using Base.Threads
using Distributed

"""
     run()
Perform the operations necessary to read in, process, and extract data from the
selected input file.

No input variables are required.

Two variables are output:
1. flags = The calculation flags from the input file.
2. coord = The molecular coordinates.

Thus, proper use of the Input.run() function would look like this:

```
flags, coord = Input.run()
```
"""
function run()
    #read in .inp and .dat files
    comm=MPI.COMM_WORLD

    if (MPI.Comm_rank(comm) == 0)
        println("-------------------------------------------------------------------------------------")
        println("                       ========================================          ")
        println("                                READING INPUT DATA FILE                  ")
        println("                       ========================================          ")
        println(" ")
    end

    directory::String = pwd()
    #println("Input file: ", directory*"/"*input_file)
    if (MPI.Comm_rank(comm) == 0)
        println(" ")
        println("Number of worker processes: ", MPI.Comm_size(comm))
        println("Number of threads per process: ", Threads.nthreads())
        println("Number of threads in total: ", MPI.Comm_size(comm)*Threads.nthreads())
    end

    flags::Flags = input_flags()

    coord::Array{Float64,2} = input_coord()
    if (MPI.Comm_rank(comm) == 0)
        println(" ")
        println("                       ========================================          ")
        println("                                       END INPUT                         ")
        println("                       ========================================          ")
    end

    return (flags,coord)
end
export run

end
