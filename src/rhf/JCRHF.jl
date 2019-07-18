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
        #== write quartet eri lists to database ==#
        eri_array::Array{Float64,1} = molecule["tei"]
        nsh::Int64 = length(basis.shells)

        eri_start::Int64 = 1
        for ish::Int64 in 1:nsh, jsh::Int64 in 1:ish
          ijsh::Int64 = index(ish,jsh)
          qnum_ij = ish*(ish-1)/2 + jsh

          ibas = basis.shells[ish].nbas
          jbas = basis.shells[jsh].nbas

          ipos = basis.shells[ish].pos
          jpos = basis.shells[jsh].pos

          for ksh::Int64 in 1:nsh, lsh::Int64 in 1:ksh
            klsh::Int64 = index(ksh,lsh)
            if (klsh > ijsh) continue end

            kbas = basis.shells[ksh].nbas
            lbas = basis.shells[lsh].nbas

            kpos = basis.shells[ksh].pos
            lpos = basis.shells[lsh].pos

            qnum_kl::Int64 = ksh*(ksh-1)/2 + lsh
            quartet_num::Int64 = qnum_ij*(qnum_ij-1)/2 + qnum_kl

            qint_ij::Int64 = ibas*(ibas-1)/2 + jbas
            qint_kl::Int64 = kbas*(kbas-1)/2 + lbas

            #eri_start::Int64 += qint_ij*(qint_ij-1)/2 + qint_kl
            #push!(eri_start_index, eri_start)

            eri_size = 0
            for μμ::Int64 in ipos:ipos+(ibas-1), νν::Int64 in jpos:jpos+(jbas-1)
              μ::Int64, ν::Int64 = μμ,νν
              if (μμ < νν) continue end

              μν::Int64 = index(μμ,νν)

              for λλ::Int64 in kpos:kpos+(kbas-1), σσ::Int64 in lpos:lpos+(lbas-1)
                λ::Int64, σ::Int64 = λλ,σσ
                if (λλ < σσ) continue end

                λσ::Int64 = index(λλ,σσ)

                if (μν < λσ)
                  two_shell::Bool = ibas == jbas
                  two_shell = two_shell || (ibas == kbas)
                  two_shell = two_shell || (ibas == lbas)
                  two_shell = two_shell || (jbas == kbas)
                  two_shell = two_shell || (jbas == lbas)
                  two_shell = two_shell || (kbas == lbas)

                  three_shell::Bool = ibas == jbas && jbas == kbas
                  three_shell = three_shell || (ibas == jbas && jbas == lbas)
                  three_shell = three_shell || (ibas == kbas && kbas == lbas)
                  three_shell = three_shell || (jbas == kbas && kbas == lbas)

                  four_shell::Bool = ibas == jbas
                  four_shell = four_shell && (jbas == kbas)
                  four_shell = four_shell && (kbas == lbas)

                  if four_shell
                      three_same::Bool = ish == jsh && jsh == ksh
                      three_same = three_same || (ish == jsh && jsh == lsh)
                      three_same = three_same || (ish == ksh && ksh == lsh)
                      three_same = three_same || (jsh == ksh && ksh == lsh)

                      four_same::Bool = ish == jsh
                      four_same = four_same && jsh == ksh
                      four_same = four_same && ksh == lsh

                      if four_same
                        if (μμ != νν && μμ != λλ && μμ != σσ &&
                          νν != λλ && νν != σσ && λλ != σσ)
                          μ,ν,λ,σ = λλ,σσ,μμ,νν
                        else
                          continue
                        end
                      elseif three_same
                        if (μμ != λλ && νν != σσ)
                          μ,ν,λ,σ = λλ,σσ,μμ,νν
                        else
                          continue
                        end
                      else
                        continue
                      end
                  elseif three_shell
                    if (μμ != λλ && νν != σσ)
                        μ,ν,λ,σ = λλ,σσ,μμ,νν
                    else
                      continue
                    end
                  elseif two_shell
                    if (μμ != λλ && νν != σσ)
                      μ,ν,λ,σ = λλ,σσ,μμ,νν
                    else
                      continue
                    end
                  end
                end

                eri_size += 1
              end
            end

            write(file, "Integrals/$quartet_num",
              eri_array[eri_start:eri_start+(eri_size-1)])

            eri_start += eri_size

            #println("$ish, $jsh, $ksh, $lsh, $quartet_num, $eri_size")
          end
        end
      end
    end
  end

  #GC.enable(false)
  scf = rhf_energy(basis, molecule, scf_flags)
  #GC.enable(true)
  #GC.gc()

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
