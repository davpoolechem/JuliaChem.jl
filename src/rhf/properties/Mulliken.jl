using JuliaChem.JCModules

function compute_mulliken_charges(mol::Molecule, basis::Basis, 
  D::Matrix{Float64}, S::Matrix{Float64})

  mulliken_charges = zeros(Float64, length(mol))
  DS_prod = D*S

  for (atom_idx, atom) in enumerate(mol)
    println("ATOM $atom_idx:")
    atom_mulliken_charge = atom.atom_id 
    println("  ATOM MULLIKEN CHARGE: $atom_mulliken_charge")
    for (shell_idx, shell) in enumerate(basis)
      println("  SHELL: $shell_idx, $(shell.atom_id)")
      if shell.atom_id == atom_idx
        atom_mulliken_charge -= DS_prod[shell_idx, shell_idx]
        println("  ATOM MULLIKEN CHARGE: $atom_mulliken_charge")
      end
      mulliken_charges[atom_idx] = atom_mulliken_charge
    end
  end     
  
  return mulliken_charges
end
