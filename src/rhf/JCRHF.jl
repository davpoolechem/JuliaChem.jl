#Base.include(@__MODULE__,"../basis/BasisStructs.jl")

"""
  module JCRHF
The module required for computation of the wave function using the *Restricted
Hartree-Fock* (RHF) method in a Self-Consistent Field (SCF) calculation. This
module will be used often, as the RHF wave function is often the zeroth-order
wave function for closed-shell systems.
"""
module JCRHF

Base.include(@__MODULE__,"RHFHelpers.jl")
Base.include(@__MODULE__,"RHFSCF.jl")
Base.include(@__MODULE__,"../../deps/src/simint.jl")

using MPI
using JSON

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
function run(mol, basis, keywords)
  comm=MPI.COMM_WORLD

  if MPI.Comm_rank(comm) == 0
      println("--------------------------------------------------------------------------------")
      println("                       ========================================                 ")
      println("                          RESTRICTED CLOSED-SHELL HARTREE-FOCK                  ")
      println("                       ========================================                 ")
      println("")
  end

  #== initialize scf flags ==#
  scf_flags = keywords["scf"]

  #== set up eris ==#
  #if MPI.Comm_rank(comm) == 0 && Threads.threadid() == 1
  if scf_flags["direct"] == true
  #  set_up_eri_database(basis)
  #else
    nshell_simint = SIMINT.allocate_shell_array(basis)
    for shell in basis.shells
      SIMINT.add_shell(shell)
    end

    SIMINT.normalize_shells()
    SIMINT.precompute_shell_pair_data()

    #for ishell::Int64 in 0:(nshell_simint-1)
    #  SIMINT.get_simint_shell_info(ishell)
    #end

    #end
  else
    println("Reading integrals from disk is not implemented yet!")
    throw()
  end

  #== actually perform scf calculation ==#
  #GC.enable(false)
  scf = rhf_energy(mol, basis, scf_flags)
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

end
