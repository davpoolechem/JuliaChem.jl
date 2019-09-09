"""
  module JCInput
The module required for reading in and processing the selected input file.
Import this module into the script when you need to process an input file
(which will be every single calculation).
"""
module JCInput

using MPI
using JSON
using Base.Threads
#using Distributed
using JLD

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
function run(args)
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
    input_string::String = read(input_file,String)
  close(input_file)

  #== initialize variables ==#
  molecule::Dict{String,Any} = Dict([])
  driver::String = ""
  model::Dict{String,Any} = Dict([])
  keywords::Dict{String,Any} = Dict([])

  #== do extraction ==#
  json_parse::Dict{String,Any} = JSON.parse(input_string)

  merge!(molecule,Dict("geometry" => json_parse["molecule"]["geometry"]))
  merge!(molecule,Dict("symbols" => json_parse["molecule"]["symbols"]))
  merge!(molecule,Dict("molecular_charge" => json_parse["molecule"]["molecular_charge"]))
  merge!(molecule,Dict("enuc" => json_parse["molecule"]["enuc"]))
  merge!(molecule,Dict("ovr" => json_parse["molecule"]["ovr"]))
  merge!(molecule,Dict("hcore" => json_parse["molecule"]["hcore"]))

  if (MPI.Comm_rank(comm) == 0) && (Threads.threadid() == 1)
    jldopen("tei_all.jld", "w") do file
      eri_array::Vector{Float64} = json_parse["molecule"]["tei"]
      write(file, "Integrals/All",eri_array)
    end
  end

  driver = json_parse["driver"]
  merge!(model,json_parse["model"])
  merge!(keywords,json_parse["keywords"])

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
