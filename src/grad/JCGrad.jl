#Base.include(@__MODULE__,"../basis/BasisStructs.jl")

"""
  module JCGrad
"""
module JCGrad

using JCModules.SIMINT
using JCModules.MolStructs
using JCModules.BasisStructs

using MPI
using JSON

function run(mol::MolStructs.Molecule, basis::BasisStructs.Basis; 
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

  nuc_grad = zeros(Float64,(length(mol.atoms),3))
  for iatom in 1:length(mol.atoms), jatom in 1:length(mol.atoms)
    if iatom != jatom
      ix = mol.atoms[iatom].atom_center[1] 
      jx = mol.atoms[jatom].atom_center[1] 

      iy = mol.atoms[iatom].atom_center[2] 
      jy = mol.atoms[jatom].atom_center[2] 

      iz = mol.atoms[iatom].atom_center[3]
      jz = mol.atoms[jatom].atom_center[3]
 
      iatm = mol.atoms[iatom].atom_id
      jatm = mol.atoms[jatom].atom_id
    
      distance = √((jx-ix)^2 + (jy-iy)^2 + (jz-iz)^2) 
    
      nuc_grad[jatom,1] -= iatm*jatm*(jx-ix)/(distance^3)
      nuc_grad[jatom,2] -= iatm*jatm*(jy-iy)/(distance^3)
      nuc_grad[jatom,3] -= iatm*jatm*(jz-iz)/(distance^3)
    end
  end 
  
  display(nuc_grad)

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
