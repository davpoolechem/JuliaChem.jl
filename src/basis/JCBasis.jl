"""
  module JCInput
The module required for reading in and processing the selected input file.
Import this module into the script when you need to process an input file
(which will be every single calculation).
"""
module JCBasis

using JCStructs

using MPI
using Base.Threads
#using Distributed
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
function run(molecule::Dict{String,Any}, model::Dict{String,Any})
  comm=MPI.COMM_WORLD

  if (MPI.Comm_rank(comm) == 0)
    println("--------------------------------------------------------------------------------")
    println("                       ========================================                 ")
    println("                                GENERATING BASIS SET                            ")
    println("                       ========================================                 ")
    println(" ")
  end

  #== initialize variables ==#
  geometry_array::Array{Float64,1} = molecule["geometry"]
  symbols::Array{String} = molecule["symbols"]
  basis::String = model["basis"]

  geometry::Array{Float64,2} = transpose(reshape(geometry_array,(3,2)))

  basis_set::Basis = Basis()
  basis_norb::Int64 = 0
  basis_nels::Int64 = 0

  #== create basis set ==#
  for atom_idx::Int64 in 1:length(symbols)
    #== initialize variables needed for shell ==#
    atom_center::Array{Float64,1} = geometry[atom_idx,:]

    #== process H shells shells ==#
    if (symbols[atom_idx] == "H")
      basis_nels += 1
      if (basis == "STO-3G")
        #== first STO3G H shell ==#
        shell_am::Int64 = 1
        shell::Shell = Shell(atom_idx, atom_center, shell_am)
        add_shell(basis_set,deepcopy(shell))
        basis_norb += 1
      end
    end
  end

  basis_set.norb = basis_norb
  basis_set.nels = basis_nels

#=
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

  =#

  if (MPI.Comm_rank(comm) == 0)
    println(" ")
    println("                       ========================================                 ")
    println("                                       END BASIS                                ")
    println("                       ========================================                 ")
  end

  return basis_set
end
export run

end
