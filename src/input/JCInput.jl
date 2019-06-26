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
#using Distributed
#using HDF5

"""
  run(args::String)
Perform the operations necessary to read in, process, and extract data from the
selected input file.

One input variable is required:
1. args = The name of the input file.

Two variables are output:
1. input_info = Information gathered from the input file.
2. basis = The basis set shells, determined from the input file.

Thus, proper use of the Input.run() function would look like this:

```
input_info, basis = Input.run(args)
```
"""
function run(args::String)
  comm=MPI.COMM_WORLD

  if (MPI.Comm_rank(comm) == 0)
    println("--------------------------------------------------------------------------------")
    println("                       ========================================                 ")
    println("                                READING INPUT DATA FILE                         ")
    println("                       ========================================                 ")
    println(" ")
  end

  #== output parallelization information ==#
  directory::String = pwd()
  #println("Input file: ", directory*"/"*input_file)
  if (MPI.Comm_rank(comm) == 0)
    println(" ")
    println("Number of worker processes: ", MPI.Comm_size(comm))
    println("Number of threads per process: ", Threads.nthreads())
    println("Number of threads in total: ",
      MPI.Comm_size(comm)*Threads.nthreads())
  end


  #== read in input file ==#
  input_file::IOStream = open(args)
    input_string::Array{String,1} = readlines(input_file)
  close(input_file)

  #== initialize variables ==#
  molecule::Dict{String,Any} = Dict([])
  driver::String = ""
  model::Dict{String,Any} = Dict([])
  keywords::Dict{String,Any} = Dict([])

  #== reformat input file and extract information from input ==#
  i::Int64 = 1
  while (i <= length(input_string))
    if (input_string[i] != "{")
      i += 1
    else
      #== do reformat ==#
      j::Int64 = i
      input_json::String = ""
      input_name::String = ""
      while (input_string[j] != "}")
        input_json *= input_string[j]
        j += 1
      end
      input_json *= "}"

      #== do extraction ==#
      json_parse::Dict{String,Any} = JSON.parse(input_json)

      merge!(molecule,json_parse["molecule"])
      driver = json_parse["driver"]
      merge!(model,json_parse["model"])
      merge!(keywords,json_parse["keywords"])

      i = j
    end
  end

  if (MPI.Comm_rank(comm) == 0)
    println(" ")
    println("                       ========================================                 ")
    println("                                       END INPUT                                ")
    println("                       ========================================                 ")
  end

  return (molecule, driver, model, keywords)
end
export run

end
