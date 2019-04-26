"""
     module JCInput
The module required for reading in and processing the selected input file.
Import this module into the script when you need to process an input file
(which will be every single calculation).
"""
module JCInput

using JCStructs

using MPI
using JSON
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
function run(args::String)
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

    if (MPI.Comm_rank(comm) == 0)
        println(" ")
        println("                       ========================================          ")
        println("                                       END INPUT                         ")
        println("                       ========================================          ")
    end

    input_file::IOStream = open(args)
        input_string::Array{String,1} = readlines(input_file)
    close(input_file)

    input_info::Dict{String,Dict{String,Any}} = Dict([])
    i::UInt32 = 1
    while (i <= length(input_string))
        if (input_string[i] != "{")
            i += 1
        else
            j::UInt32 = i
            json_subsection::String = ""
            input_name::String = ""
            while (input_string[j] != "}")
                json_subsection *= input_string[j]
                if (occursin(r"\\\"Input\\\"\:(.*)",input_string[j]))
                    input_name = match(r"\\\"Input\\\"\:(.*)",input_string[j])[1]
                    input_name = input_name[2:end-2]
                end
                j += 1
            end
            json_subsection *= "}"
            json_parse::Dict{String,Any} = JSON.parse(json_subsection)
            merge!(input_info,Dict([(input_name,json_parse)]))

            i = j
        end
    end

    shell_am = input_info["Basis Flags"]["shells"]
    basis::Basis = Basis()
    for i in 1:length(shell_am)
        shell::Shell = Shell(UInt32(shell_am[i]))
        add_shell(basis,deepcopy(shell))
    end

    return (input_info,basis)
end
export run

end
