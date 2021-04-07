using JuliaChem.JCModules

using HDF5

function compute_formation_energy(mol::Molecule, basis::Basis,
  total_energy::Float64)

  #== generate list of unique atoms in system ==#
  atom_list = Dict{String, Float64}([])
  h5open(joinpath(@__DIR__, "../../../records/eatom.h5"),"r") do eatom               
    for atom in mol
      if !haskey(atom_list, atom.symbol)
        atom_list[atom.symbol] = read(eatom[atom.symbol][basis.model])       
      end
    end      
  end     
 
  #== compute energy sum of constituent atoms ==#
  constituent_energy_sum = 0.0
  for atom in mol
    constituent_energy_sum += atom_list[atom.symbol]
  end 

  #== compute formation energy ==#
  E_form = total_energy - constituent_energy_sum

  #== we are done ==#
  return E_form 
end
