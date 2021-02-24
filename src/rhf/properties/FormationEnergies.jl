using JuliaChem.JCModules

function compute_formation_energy(mol::Molecule, basis::Basis,
  total_energy::Float64)

  #== generate list of unique atoms in system ==#
  h5open(joinpath(@__DIR__, "../../../records/eatom.h5"),"r") do eatom               
    atom_list = Dict{String, Float64}([])
    
    for atom in mol
      if !haskey(atom_list, atom)
        atom_list[atom] = eatom[atom][basis]       
      end
    end      
  end     
 
  #== compute energy sum of constituent atoms ==#
  constituent_energy_sum = 0.0
  for atom in mol
    constituent_energy_sum += atom_list[atom]
  end 

  #== compute formation energy ==#
  E_form = total_energy - constituent_energy_sum

  #== we are done ==#
  return E_form 
end
