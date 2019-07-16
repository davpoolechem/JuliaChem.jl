"""
  module JCRHF
The module required for computation of the wave function using the *Restricted
Hartree-Fock* (RHF) method in a Self-Consistent Field (SCF) calculation. This
module will be used often, as the RHF wave function is often the zeroth-order
wave function for closed-shell systems.
"""
module JCRHF

Base.include(@__MODULE__,"RHFSCF.jl")

using BasisStructs

using MPI
using JSON
using HDF5

"""
  run(input_info::Dict{String,Dict{String,Any}}, basis::Basis)

Execute the JuliaChem RHF algorithm.

One input variable is required:
1. input_info = Information gathered from the input file.
2. basis = The basis set shells, determined from the input file.

One variable is output:
1. scf = Data saved from the SCF calculation.

Thus, proper use of the RHF.run() function would look like this:

```
scf = RHF.run(input_info, basis)
```
"""
function run(basis::Basis, molecule::Dict{String,Any},
  keywords::Dict{String,Any})

  comm=MPI.COMM_WORLD

  if (MPI.Comm_rank(comm) == 0)
      println("--------------------------------------------------------------------------------")
      println("                       ========================================                 ")
      println("                          RESTRICTED CLOSED-SHELL HARTREE-FOCK                  ")
      println("                       ========================================                 ")
      println("")
  end

  #== initialize scf flags ==#
  scf_flags::Dict{String,Any} = keywords["scf"]

  #== set up eri database if not doing direct ==#
  if (scf_flags["direct"] == false)
    hdf5name = "tei"
    hdf5name *= ".h5"
    if ((MPI.Comm_rank(comm) == 0) && (Threads.threadid() == 1))
      h5open(hdf5name, "w") do file
        eri_array::Array{Float64,1} = molecule["tei"]
        write(file, "tei", eri_array)

        nsh::Int64 = length(basis.shells)

        eri_start_index::Array{Int64,1} = [ ]
        eri_index_size::Int64 = 0
        for ish::Int64 in 1:nsh, jsh::Int64 in 1:ish
          ijsh::Int64 = index(ish,jsh)
          #qnum_ij = ish*(ish-1)/2 + jsh

          ibas = basis.shells[ish].nbas
          jbas = basis.shells[jsh].nbas

          for ksh::Int64 in 1:nsh, lsh::Int64 in 1:ksh
            klsh::Int64 = index(ksh,lsh)
            if (klsh > ijsh) continue end

            kbas = basis.shells[ksh].nbas
            lbas = basis.shells[lsh].nbas

            #qnum_kl::Int64 = ksh*(ksh-1)/2 + lsh
            #quartet_num::Int64 = qnum_ij*(qnum_ij-1)/2 + qnum_kl

            qint_ij::Int64 = ibas*(ibas-1)/2 + jbas
            qint_kl::Int64 = kbas*(kbas-1)/2 + lbas

            eri_index_size::Int64 += qint_ij*(qint_ij-1)/2 + qint_kl
            push!(eri_start_index, eri_index_size)
          end
        end

        write(file, "start", eri_start_index)
      end
    end
  end

  #GC.enable(false)
  scf = rhf_energy(basis, molecule, scf_flags)
  #GC.enable(true)
  #GC.gc()
#=
  calculation_name::String = ctrl_flags.NAME
  json_output = open("output.json","w")
    output_name = Dict([("Calculation",calculation_name)])
    output_fock = Dict([("Structure","Fock"),("Data",scf.Fock)])
    output_density = Dict([("Structure","Density"),("Data",scf.Density)])
    output_coeff = Dict([("Structure","Coeff"),("Data",scf.Coeff)])
    if (MPI.Comm_rank(comm) == 0)
      write(json_output,JSON.json(output_name))
      write(json_output,JSON.json(output_fock))
      write(json_output,JSON.json(output_density))
      write(json_output,JSON.json(output_coeff))
    end
  close(json_output)
=#
  if (MPI.Comm_rank(comm) == 0)
    println("                       ========================================                 ")
    println("                             END RESTRICTED CLOSED-SHELL                 ")
    println("                                     HARTREE-FOCK                        ")
    println("                       ========================================                 ")
  end

  return scf
end
export run

#=
function run(flags::RHF_Flags, restart::RHFRestartData)
    comm=MPI.COMM_WORLD

    if (MPI.Comm_rank(comm) == 0)
        println("--------------------------------------------------------------------------------------")
        println("                       ========================================          ")
        println("                         RESTRICTED CLOSED-SHELL HARTREE-FOCK            ")
        println("                       ========================================          ")
        println("")
    end

    #GC.enable(false)
    scf = rhf_energy(flags,restart)
    #GC.enable(true)
    #GC.gc()

    if (MPI.Comm_rank(comm) == 0)
        println("                       ========================================          ")
        println("                             END RESTRICTED CLOSED-SHELL                 ")
        println("                                     HARTREE-FOCK                        ")
        println("                       ========================================          ")
    end

    output_fock = Dict([("Structure","Fock"),("Data",scf.Fock)])
    output_density = Dict([("Structure","Density"),("Data",scf.Density)])
    output_coeff = Dict([("Structure","Coeff"),("Data",scf.Coeff)])
    if (MPI.Comm_rank(comm) == 0)
        json_output = open("test.json","w")
            write(json_output,JSON.json(output_fock))
            write(json_output,JSON.json(output_density))
            write(json_output,JSON.json(output_coeff))
        close(json_output)
    end

    return scf
end
export run
=#

end
