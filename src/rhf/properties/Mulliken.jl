using JuliaChem.JCModules

function compute_population_analysis(mol::Molecule, basis::Basis, 
  D::Matrix{Float64}, S::Matrix{Float64})

  #== compute mulliken/lowdin population per ao ==#
  ao_mulliken_pop = D.*S

  #S_one_half = S^0.5
  #ao_lowdin_pop = S_one_half.*D.*S_one_half

  #== compute mulliken/lowdin population per atom ==#
  atom_mulliken_overlap_pop = zeros(Float64, (length(mol), length(mol))) 
  #atom_lowdin_overlap_pop = zeros(Float64, (length(mol), length(mol))) 
  for iatom in 1:length(mol), jatom in 1:length(mol)
    iatom_aos = [] 
    jatom_aos = [] 
   
    for shell in basis
      if shell.atom_id == iatom 
        push!(iatom_aos, shell.pos + shell.nbas - 1) 
      end
      if shell.atom_id == jatom 
        push!(jatom_aos, shell.pos + shell.nbas - 1) 
      end
    end

    iatom_start = iatom_aos[1]
    iatom_end = iatom_aos[end] 

    jatom_start = jatom_aos[1]
    jatom_end = jatom_aos[end] 

    for iao in iatom_start:iatom_end, jao in jatom_start:jatom_end
      i = max(iao, jao)
      j = min(iao, jao)
      atom_mulliken_overlap_pop[iatom, jatom] += ao_mulliken_pop[i,j] 
      #atom_lowdin_overlap_pop[iatom, jatom] += ao_lowdin_pop[i,j] 
    end
  end 

  atom_mulliken_pop = zeros(Float64, length(mol)) 
  #atom_lowdin_pop = zeros(Float64, length(mol)) 
  for iatom in 1:length(mol)
    @views atom_mulliken_pop[iatom] = sum(atom_mulliken_overlap_pop[:,iatom])
    #@views atom_lowdin_pop[iatom] = sum(atom_lowdin_overlap_pop[:,iatom])
  end

  #for (atom_idx, atom) in enumerate(mol)
  #  atom_mulliken_charge = atom.atom_id 
  #  for (shell_idx, shell) in enumerate(basis)
  #    if shell.atom_id == atom_idx
  #      atom_mulliken_ -= DS_prod[shell_idx, shell_idx]
  #    end
  #    mulliken_charges[atom_idx] = atom_mulliken_charge
  #  end
  #end     
  
  return atom_mulliken_pop #, ao_lowdin_pop
end
