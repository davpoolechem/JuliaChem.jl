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
using HDF5

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
  #read in .inp and .dat files
  comm=MPI.COMM_WORLD

  if (MPI.Comm_rank(comm) == 0)
    println("--------------------------------------------------------------------------------")
    println("                       ========================================                 ")
    println("                                GENERATING BASIS SET                            ")
    println("                       ========================================                 ")
    println(" ")
  end


  if (MPI.Comm_rank(comm) == 0)
    println(" ")
    println("                       ========================================                 ")
    println("                                       END BASIS                                ")
    println("                       ========================================                 ")
  end

  shell_am = input_info["Basis Flags"]["shells"]
  basis::Basis = Basis()
  for i in 1:length(shell_am)
    shell::Shell = Shell(UInt32(shell_am[i]))
    add_shell(basis,deepcopy(shell))
  end

  #== set up eri database ==#
  norb = input_info["Basis Flags"]["norb"]

  hdf5name = input_info["Control Flags"]["name"]
  hdf5name *= ".h5"
  if ((MPI.Comm_rank(comm) == 0) && (Threads.threadid() == 1))
    h5open(hdf5name, "w") do file
      eri_array::Array{Float64,1} = input_info["Two-Electron"]["tei"]
      write(file, "tei", eri_array)
    end
  end
  MPI.Barrier

  return basis
end
export run

end
