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

  #== compute nuclear gradient ==#
  nuc_grad = compute_nuc_grad(mol) 
  println("NUC GRAD")
  display(nuc_grad); println()

  #== compute one-electron gradient ==#
  S_grad = zeros(Float64, (basis.norb, basis.norb))
  compute_overlap_grad(S_grad, basis, jeri_oei_grad_engine) 
  println("S GRAD")
  display(S_grad); println()
 
  if MPI.Comm_rank(comm) == 0 && output == "verbose"
    println("                       ========================================                 ")
    println("                              END RESTRICTED CLOSED-SHELL                       ")
    println("                                 HARTREE-FOCK GRADIENT                          ")
    println("                       ========================================                 ")
    println("--------------------------------------------------------------------------------")
  end

  return nuc_grad 
end
export run

end
