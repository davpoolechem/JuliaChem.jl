#Base.include(@__MODULE__,"../basis/BasisStructs.jl")

"""
  module JCRHF
The module required for computation of the wave function using the *Restricted
Hartree-Fock* (RHF) method in a Self-Consistent Field (SCF) calculation. This
module will be used often, as the RHF wave function is often the zeroth-order
wave function for closed-shell systems.
"""
module JCRHF

using JCModules.BasisStructs

Base.include(@__MODULE__,"RHFSCF.jl")

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
function run(basis::BasisStructs.Basis, molecule::Dict{String,Any},
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
    if ((MPI.Comm_rank(comm) == 0) && (Threads.threadid() == 1))
      h5open("tei_batch.h5", "w") do file
        #== write quartet eri lists to database ==#
        eri_array::Array{Float64,1} = []
        h5open("tei_all.h5", "r") do tei
          eri_array = read(tei,"Integrals/All")
        end

        nsh::Int64 = length(basis.shells)

        eri_array_batch::Array{Float64,1} = [ ]
        eri_array_starts::Array{Float64,1} = [ ]
        eri_array_sizes::Array{Float64,1} = [ ]

        eri_start::Int64 = 1
        quartet_batch_num_old::Int64 = 1

        for ish::Int64 in 1:nsh, jsh::Int64 in 1:ish
          ijsh::Int64 = index(ish,jsh)
          qnum_ij = ish*(ish-1)/2 + jsh

          ibas::Int64 = basis.shells[ish].nbas
          jbas::Int64 = basis.shells[jsh].nbas

          ipos::Int64 = basis.shells[ish].pos
          jpos::Int64 = basis.shells[jsh].pos

          for ksh::Int64 in 1:nsh, lsh::Int64 in 1:ksh
            klsh::Int64 = index(ksh,lsh)
            if (klsh > ijsh) continue end

            kbas::Int64 = basis.shells[ksh].nbas
            lbas::Int64 = basis.shells[lsh].nbas

            kpos::Int64 = basis.shells[ksh].pos
            lpos::Int64 = basis.shells[lsh].pos

            qnum_kl::Int64 = ksh*(ksh-1)/2 + lsh
            quartet_num::Int64 = qnum_ij*(qnum_ij-1)/2 + qnum_kl

            qint_ij::Int64 = ibas*(ibas-1)/2 + jbas
            qint_kl::Int64 = kbas*(kbas-1)/2 + lbas

            #eri_start::Int64 += qint_ij*(qint_ij-1)/2 + qint_kl
            #push!(eri_start_index, eri_start)

            eri_size::Int64 = 0
            println("$ish, $jsh, $ksh, $lsh")
            for μμ::Int64 in ipos:ipos+(ibas-1), νν::Int64 in jpos:jpos+(jbas-1)
              μ::Int64, ν::Int64 = μμ,νν
              if (μμ < νν) continue end

              μν::Int64 = index(μμ,νν)

              for λλ::Int64 in kpos:kpos+(kbas-1), σσ::Int64 in lpos:lpos+(lbas-1)
                λ::Int64, σ::Int64 = λλ,σσ
                if (λλ < σσ) continue end

                λσ::Int64 = index(λλ,σσ)

                #print("$μμ, $νν, $λλ, $σσ;") 
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
                        #print("\n")
                        continue
                      elseif three_same
                        if (μμ != λλ && νν != σσ)
                          μ,ν,λ,σ = λλ,σσ,μμ,νν
                        else
                          #print("\n")
                          continue
                        end
                      else
                          #print("\n")
                        continue
                      end
                  elseif three_shell
                    if (μμ != λλ && νν != σσ)
                        μ,ν,λ,σ = λλ,σσ,μμ,νν
                    else
                          #print("\n")
                      continue
                    end
                  elseif two_shell
                    if (ish == ksh && jsh == lsh &&
                        μμ != νν && μμ != λλ && μμ != σσ &&
                        νν != λλ && νν != σσ &&
                        λλ != σσ )
                         #print("\n")
                      continue
                    elseif (μμ != λλ && νν != σσ)
                      μ,ν,λ,σ = λλ,σσ,μμ,νν
                    else
                          #print("\n")
                      continue
                    end
                  end
                end
                #println(" $μ, $ν, $λ, $σ")
                eri_size += 1
              end
            end

            quartets_per_batch::Int64 = 2500
            quartet_batch_num::Int64 = Int64(floor(quartet_num/
              quartets_per_batch)) + 1

            display(eri_array[eri_start:eri_start+(eri_size-1)])
            println(" ")
            if quartet_batch_num != quartet_batch_num_old
              #== write arrays to disk ==#
              write(file, "Integrals/$quartet_batch_num_old",
                eri_array_batch)
              write(file, "Starts/$quartet_batch_num_old",
                eri_array_starts)
              write(file, "Sizes/$quartet_batch_num_old",
                eri_array_sizes)

              #== reset variables as needed ==#
              eri_array_batch = [ ]
              eri_array_starts = [ ]
              eri_array_sizes = [ ]

              quartet_batch_num_old = quartet_batch_num
            else
              eri_array_batch = [ eri_array_batch;
              eri_array[eri_start:eri_start+(eri_size-1)]]

              eri_start_readin::Int64 = eri_start - quartets_per_batch*
                (quartet_batch_num-1)
              push!(eri_array_starts,eri_start_readin)

              push!(eri_array_sizes,eri_size)
            end

            eri_start += eri_size

            #println("$ish, $jsh, $ksh, $lsh, $quartet_num, $eri_size")
          end
        end

        write(file, "Integrals/$quartet_batch_num_old",
          eri_array_batch)
        write(file, "Starts/$quartet_batch_num_old",
          eri_array_starts)
        write(file, "Sizes/$quartet_batch_num_old",
          eri_array_sizes)
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
