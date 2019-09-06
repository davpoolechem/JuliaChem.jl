#Base.include(@__MODULE__,"../basis/BasisStructs.jl")

"""
  module JCInput
The module required for reading in and processing the selected input file.
Import this module into the script when you need to process an input file
(which will be every single calculation).
"""
module JCBasis

using JCModules.BasisStructs

using MPI
using Base.Threads
#using Distributed
using HDF5

Base.include(@__MODULE__, "BasisHelpers.jl")

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
function run(molecule, model)
  comm=MPI.COMM_WORLD

  if (MPI.Comm_rank(comm) == 0)
    println("--------------------------------------------------------------------------------")
    println("                       ========================================                 ")
    println("                                GENERATING BASIS SET                            ")
    println("                       ========================================                 ")
    println(" ")
  end

  #== initialize variables ==#
  geometry_array::Vector{Float64} = molecule["geometry"]
  symbols::Vector{String} = molecule["symbols"]
  basis::String = model["basis"]
  charge::Int64 = molecule["molecular_charge"]

  num_atoms::Int64 = length(geometry_array)/3
  geometry_array_t::Matrix{Float64} = reshape(geometry_array,(3,num_atoms))
  geometry::Matrix{Float64} = transpose(geometry_array_t)

  basis_set::Basis = Basis(basis, charge)
  atomic_number_mapping::Dict{String,Int64} = create_atomic_number_mapping()
  shell_am_mapping::Dict{String,Int64} = create_shell_am_mapping()

  println("----------------------------------------          ")
  println("        Basis Set Information...                  ")
  println("----------------------------------------          ")

  #== create basis set ==#
  hdf5name::String = "bsed"
  hdf5name *= ".h5"
  h5open(hdf5name,"r") do bsed
    for atom_idx::Int64 in 1:length(symbols)
      #== initialize variables needed for shell ==#
      atom_center::Vector{Float64} = geometry[atom_idx,:]

      symbol::String = symbols[atom_idx]
      atomic_number::Int64 = atomic_number_mapping[symbol]

      basis_set.nels += atomic_number

      #== read in basis set values==#
      shells::Dict{String,Any} = read(
        bsed["$symbol/$basis"])

      #== process basis set values into shell objects ==#
      println("ATOM $symbol:")
      for shell_num::Int64 in 1:length(shells)
        new_shell_dict::Dict{String,Any} = shells["$shell_num"]
        println("Shell #$shell_num:")
        display(new_shell_dict)
        println("")

        new_shell_am::Int64 = shell_am_mapping[new_shell_dict["Shell Type"]]
        new_shell_exp::Vector{Float64} = new_shell_dict["Exponents"]
        new_shell_coeff::Array{Float64} = new_shell_dict["Coefficients"]

        new_shell_nprim::Int64 = size(new_shell_exp)[1]
        new_shell_coeff_array::Vector{Float64} = reshape(new_shell_coeff,
          (length(new_shell_coeff),))       
  
        new_shell = Shell(atom_idx, new_shell_exp, new_shell_coeff_array,
          atom_center, new_shell_am, size(new_shell_exp)[1])
        add_shell(basis_set,deepcopy(new_shell))

        basis_set.norb += new_shell.nbas
      end
      println(" ")
      println(" ")

    end
  end

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
