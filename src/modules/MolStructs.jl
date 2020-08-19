struct Atom 
  atom_id::Int64
  symbol::String

  atom_center::SVector{3,Float64}
end
export Atom 

struct Molecule
  atoms::Vector{Atom}
end
export Molecule 

function Base.getindex(mol::Molecule, index)
  return mol.atoms[index]
end

