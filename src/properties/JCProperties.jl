#Base.include(@__MODULE__,"../basis/BasisStructs.jl")

"""
  module JCProperties
"""
module JCProperties 

include("PropHelpers.jl")

using JuliaChem.JCModules
using JuliaChem.JERI

using JSON
using MPI
using Printf

function run(mol::Molecule, basis::Basis, rhf_energy, 
  keywords; output="none")
  
  comm=MPI.COMM_WORLD

  if MPI.Comm_rank(comm) == 0 && output == "verbose"
      println("--------------------------------------------------------------------------------")
      println("                       ========================================                 ")
      println("                                RESTRICTED CLOSED-SHELL                         ")
      println("                                HARTREE-FOCK PROPERTIES                         ")
      println("                       ========================================                 ")
      println("")
  end

  #== initial setup ==#
  jeri_prop_engine = JERI.PropEngine(mol.mol_cxx, 
    basis.basis_cxx) 

  P = rhf_energy["Density"] 

  #== compute dipole moment ==#
  if MPI.Comm_rank(comm) == 0 && output == "verbose"
    println("----------------------------------------          ")
    println("     Computing multiple moments...                ")
    println("----------------------------------------          ")
    println(" ")
    println("Dipole:       X           Y           Z         Tot. (D)        ") 
  end  
  
  dipole = compute_dipole(mol, basis, P, jeri_prop_engine)
  dipole_moment = sqrt(dipole[1]^2 + dipole[2]^2 + dipole[3]^2)
  
  @printf("          %.6f   %.6f    %.6f    %.6f     \n", 
    dipole[1], dipole[2], dipole[3], dipole_moment)

  if MPI.Comm_rank(comm) == 0 && output == "verbose"
    println("                       ========================================                 ")
    println("                              END RESTRICTED CLOSED-SHELL                       ")
    println("                                HARTREE-FOCK PROPERTIES                         ")
    println("                       ========================================                 ")
    println("--------------------------------------------------------------------------------")
  end
end
export run

end
