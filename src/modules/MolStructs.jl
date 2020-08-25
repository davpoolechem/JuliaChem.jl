#== structure representing an Atom ==#
struct Atom 
  atom_id::Int64
  symbol::String

  atom_center::SVector{3,Float64}
end
export Atom 

#== structure representing a molecule ==#
struct Molecule
  atoms::Vector{Atom}
  mol_cxx
end
export Molecule 

#== functions for Molecule ==#
function Base.getindex(mol::Molecule, index)
  return mol.atoms[index]
end

function Base.length(mol::Molecule)
  return length(mol.atoms)
end

function Base.iterate(mol::Molecule)
  return iterate(mol.atoms)
end

function Base.iterate(mol::Molecule, state)
  return iterate(mol.atoms, state)
end

function Base.push!(mol::Molecule, atom::Atom)
  return push!(mol.atoms, atom)
end
