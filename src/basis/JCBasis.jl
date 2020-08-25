#Base.include(@__MODULE__,"../basis/BasisStructs.jl")
"""
  module JCInput
The module required for reading in and processing the selected input file.
Import this module into the script when you need to process an input file
(which will be every single calculation).
"""
module JCBasis

using JuliaChem.JCModules
using JuliaChem.JERI

using CxxWrap
using MPI
using Base.Threads
using HDF5
using PrettyTables

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
function run(molecule, model; output="none")
  comm=MPI.COMM_WORLD

  if MPI.Comm_rank(comm) == 0 && output == "verbose"
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

  atomic_number_mapping::Dict{String,Int64} = create_atomic_number_mapping()
  shell_am_mapping::Dict{String,Int64} = create_shell_am_mapping()

  mol = Molecule([], StdVector{JERI.Atom}())

  if MPI.Comm_rank(comm) == 0 && output == "verbose"
    println("----------------------------------------          ")
    println("        Basis Set Information...                  ")
    println("----------------------------------------          ")
  end

  basis_set_shells = Vector{JCModules.Shell}([])
  shells_cxx = StdVector([ StdVector{JERI.Shell}() for i in 1:55 ]) 
  shells_cxx_added = [ false for i in 1:55 ]

  basis_set_nels = -charge 
  basis_set_norb = 0
  pos = 1

  #== create basis set ==#
  h5open(joinpath(@__DIR__, "../../records/bsed.h5"),"r") do bsed
    shell_id = 1
    for atom_idx::Int64 in 1:length(symbols)
      #== initialize variables needed for shell ==#
      atom_center::Vector{Float64} = geometry[atom_idx,:]
      atom_center[:] .*= 1.0/0.52917724924 #switch from angs to bohr

      symbol::String = symbols[atom_idx]
      atomic_number::Int64 = atomic_number_mapping[symbol]

      #== create atom objects ==#
      push!(mol, JCModules.Atom(atomic_number, symbol, atom_center))
      push!(mol.mol_cxx, JERI.create_atom(atomic_number, atom_center)) 
       
      basis_set_nels += atomic_number

      #== read in basis set values==#
      shells::Dict{String,Any} = read(
        bsed["$symbol/$basis"])

      #== process basis set values into shell objects ==#
      if MPI.Comm_rank(comm) == 0 && output == "verbose"
        println("ATOM #$atom_idx ($symbol):") 
      end
      
      for shell_num::Int64 in 1:length(shells)
        new_shell_dict::Dict{String,Any} = shells["$shell_num"]

        new_shell_am::Int64 = shell_am_mapping[new_shell_dict["Shell Type"]]
        new_shell_exp::Vector{Float64} = new_shell_dict["Exponents"]
        new_shell_coeff::Array{Float64} = new_shell_dict["Coefficients"]

        #== if L shell, divide up ==# 
        if new_shell_am == -1
          #== s component ==#
          if MPI.Comm_rank(comm) == 0 && output == "verbose"
            println("L (s)")
            pretty_table(hcat(collect(1:length(new_shell_exp)),new_shell_exp, 
              new_shell_coeff[:,1]), 
              vcat( [ "Primitive" "Exponent" "Contraction Coefficient" ] ),
              formatters = ft_printf("%5.6f", [2,3]) )
          end 

          new_shell_nprim = size(new_shell_exp)[1]

          new_shell = JCModules.Shell(shell_id, atom_idx, atomic_number, 
            new_shell_exp, new_shell_coeff[:,1],
            atom_center, 1, size(new_shell_exp)[1], pos, true)
          push!(basis_set_shells, new_shell)
          #display(new_shell_coeff[:,1])
          if !shells_cxx_added[atomic_number+1]
            push!(shells_cxx[atomic_number+1], JERI.create_shell(0, 
              new_shell_exp, new_shell_coeff[:,1], atom_center))
          end

          basis_set_norb += 1 
          shell_id += 1
          pos += new_shell.nbas

          #== p component ==#
          if MPI.Comm_rank(comm) == 0 && output == "verbose"
            println("L (p)")
            pretty_table(hcat(collect(1:length(new_shell_exp)),new_shell_exp, 
              new_shell_coeff[:,2]), 
              vcat( [ "Primitive" "Exponent" "Contraction Coefficient" ] ),
              formatters = ft_printf("%5.6f", [2,3]) )
          end 

          new_shell_nprim = size(new_shell_exp)[1]

          new_shell = JCModules.Shell(shell_id, atom_idx, atomic_number,
            new_shell_exp, new_shell_coeff[:,2],
            atom_center, 2, size(new_shell_exp)[1], pos, true)
          push!(basis_set_shells,new_shell)
          #display(new_shell_coeff[:,2])
          if !shells_cxx_added[atomic_number+1]
            push!(shells_cxx[atomic_number+1], JERI.create_shell(1, 
              new_shell_exp, new_shell_coeff[:,2], atom_center))
          end

          basis_set_norb += 3 
          shell_id += 1
          pos += new_shell.nbas
        #== otherwise accept shell as is ==#
        else 
          if MPI.Comm_rank(comm) == 0 && output == "verbose"
            println(new_shell_dict["Shell Type"])
            pretty_table(hcat(collect(1:length(new_shell_exp)),new_shell_exp, 
              new_shell_coeff), 
              vcat( [ "Primitive" "Exponent" "Contraction Coefficient" ] ),
              formatters = ft_printf("%5.6f", [2,3]) )
          end 

          new_shell_nprim = size(new_shell_exp)[1]
          new_shell_coeff_array = reshape(new_shell_coeff,
            (length(new_shell_coeff),))       

          new_shell = JCModules.Shell(shell_id, atom_idx, atomic_number,
            new_shell_exp, deepcopy(new_shell_coeff_array),
            atom_center, new_shell_am, size(new_shell_exp)[1], pos, true)
          push!(basis_set_shells,new_shell)
          if !shells_cxx_added[atomic_number+1]
            push!(shells_cxx[atomic_number+1], JERI.create_shell(new_shell_am-1, 
              new_shell_exp, new_shell_coeff_array, atom_center))
          end 

          basis_set_norb += new_shell.nbas
          shell_id += 1
          pos += new_shell.nbas
        end

        if MPI.Comm_rank(comm) == 0 && output == "verbose"
          println(" ")
        end
      end
      
      shells_cxx_added[atomic_number+1] = true 
      #display(shells_cxx)

      if MPI.Comm_rank(comm) == 0 && output == "verbose"
        println(" ")
        println(" ")
      end
    end
  end

  #sort!(basis_set_shells, by = x->((x.nbas*x.nprim),x.am))
  #sort!(basis_set_shells, by = x->(x.atomic_number,x.atom_id))
  
  basis_set::Basis = Basis(basis_set_shells, shells_cxx, basis, 
    basis_set_norb, basis_set_nels)                                       

  #== set up shell pair ordering ==#
  #for ish in 1:length(basis_set.shells), jsh in 1:ish
  #  push!(basis_set.shpair_ordering, ShPair(basis_set.shells[ish], 
  #    basis_set.shells[jsh])) 
  #end
  
  #sort!(basis_set.shpair_ordering, by = x->((x.nbas2*x.nprim2),x.am2,unsafe_string(x.class)))
 
  #for ish in 1:length(basis_set.shells), jsh in 1:ish 
  #  idx = ceil(Int64, ish*(ish-1)/2) + jsh
  #  push!(basis_set.shpair_ordering,(shellpairs[idx].sh_a.shell_id, 
  #    shellpairs[idx].sh_b.shell_id))
  #end

  #for shellpair in shellpairs 
  #  basis_set.shpair_ordering = vcat(basis_set.shpair_ordering, [shellpair.sh_a.shell_id shellpair.sh_b.shell_id])
  #end
  
  #delete first row, as it is simply zeroes
  #basis_set.shpair_ordering = basis_set.shpair_ordering[setdiff(1:end, 1),:]

  if MPI.Comm_rank(comm) == 0 && output == "verbose"
    println(" ")
    println("                       ========================================                 ")
    println("                                       END BASIS                                ")
    println("                       ========================================                 ")
  end

  return mol, basis_set
end
export run

end
