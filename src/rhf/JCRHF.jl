#Base.include(@__MODULE__,"../basis/BasisStructs.jl")

"""
  module JCRHF
The module required for computation of the wave function using the *Restricted
Hartree-Fock* (RHF) method in a Self-Consistent Field (SCF) calculation. This
module will be used often, as the RHF wave function is often the zeroth-order
wave function for closed-shell systems.
"""
module JCRHF

using JCModules.SIMINT
using MPI
using JSON

#                        | s | p |    d     |    f    |
#                        ------------------------------ 
const axial_norm_fact = [ 1.0 1.0    1.0        1.0   ; 
                          0.0 1.0 sqrt(3.0)  sqrt(5.0);
                          0.0 1.0 sqrt(3.0)  sqrt(5.0);
                          0.0 0.0    1.0     sqrt(5.0);
                          0.0 0.0 sqrt(3.0) sqrt(15.0);
                          0.0 0.0    1.0     sqrt(5.0);
                          0.0 0.0    0.0        1.0   ;
                          0.0 0.0    0.0     sqrt(5.0);
                          0.0 0.0    0.0     sqrt(5.0);
                          0.0 0.0    0.0        1.0
                        ]

Base.include(@__MODULE__,"RHFHelpers.jl")
Base.include(@__MODULE__,"RHFSCF.jl")


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
function run(mol::MolStructs.Molecule, basis::BasisStructs.Basis, 
  scf_flags; output="none")
  
  comm=MPI.COMM_WORLD

  if MPI.Comm_rank(comm) == 0 && output == "verbose"
    println("--------------------------------------------------------------------------------")
    println("                       ========================================                 ")
    println("                                RESTRICTED CLOSED-SHELL                         ")
    println("                                  HARTREE-FOCK ENERGY                           ")
    println("                       ========================================                 ")
    println("")
  end

  #== set up eris ==#
  if scf_flags["direct"] == true
    nshell_simint = SIMINT.allocate_shell_array(basis)
    for shell in basis.shells
      SIMINT.add_shell(shell)
    end

    SIMINT.precompute_shell_pair_data()

    #for ishell::Int64 in 0:(nshell_simint-1)
    #  SIMINT.get_simint_shell_info(ishell)
    #end
  else
    println("Reading integrals from disk is not implemented yet!")
    throw()
  end

  #== actually perform scf calculation ==#
  rhfenergy = rhf_energy(mol, basis, scf_flags; output=output)

  if MPI.Comm_rank(comm) == 0 && output == "verbose"
    println("                       ========================================                 ")
    println("                              END RESTRICTED CLOSED-SHELL                       ")
    println("                                  HARTREE-FOCK ENERGY                           ")
    println("                       ========================================                 ")
    println("--------------------------------------------------------------------------------")
  end

  return rhfenergy 
end
export run

end
