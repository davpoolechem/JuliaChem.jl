using BasisStructs

#== add H shells to basis set ==#
function process_H_shell(basis_set::Basis, atom_idx::Int64,
  atom_center::Array{Float64,1})

  basis_set.nels += 1

  #== STO-3G basis ==#
  if (basis_set.model == "STO-3G")
    #== first STO-3G H shell ==#
    shell_am::Int64 = 1
    shell::Shell = Shell(atom_idx, atom_center, shell_am)
    add_shell(basis_set,deepcopy(shell))
    basis_set.norb += 1
  end
end
export process_H_shell
