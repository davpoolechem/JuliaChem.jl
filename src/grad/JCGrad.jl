#Base.include(@__MODULE__,"../basis/BasisStructs.jl")

"""
  module JCGrad
"""
module JCGrad

include("GradHelpers.jl")

using JuliaChem.JCModules
using JuliaChem.JERI

using MPI
using JSON

function run(mol::Molecule, basis::Basis, rhf_energy; 
  output="none")
  
  comm=MPI.COMM_WORLD

  if MPI.Comm_rank(comm) == 0 && output == "verbose"
      println("--------------------------------------------------------------------------------")
      println("                       ========================================                 ")
      println("                                RESTRICTED CLOSED-SHELL                         ")
      println("                                 HARTREE-FOCK GRADIENT                          ")
      println("                       ========================================                 ")
      println("")
  end

  #== initial setup ==#
  jeri_oei_grad_engine = JERI.OEIEngine(mol.mol_cxx, 
    basis.basis_cxx, 1) 
  W = rhf_energy["Energy-Weighted Density"] 

  #== compute nuclear gradient contribution==#
  nuclear_gradient = compute_nuc_grad(mol) 
  println("NUCLEAR REPULSION: ")
  display(nuclear_gradient); println()

  #== compute overlap gradient contribution ==#
  overlap_gradient = compute_overlap_grad(mol, basis, W, jeri_oei_grad_engine) 
  println("OVERLAP GRADIENT: ")
  display(overlap_gradient); println()
  
  total_gradient = nuclear_gradient .+ overlap_gradient
  println("TOTAL GRADIENT: ")
  display(total_gradient); println()
 
  if MPI.Comm_rank(comm) == 0 && output == "verbose"
    println("                       ========================================                 ")
    println("                              END RESTRICTED CLOSED-SHELL                       ")
    println("                                 HARTREE-FOCK GRADIENT                          ")
    println("                       ========================================                 ")
    println("--------------------------------------------------------------------------------")
  end

  return nuclear_gradient 
end
export run

end
