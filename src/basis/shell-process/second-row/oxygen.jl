#== add O shells to basis set ==#
function process_O_shell(basis_set::Basis, atom_idx::Int64,
  atom_center::Array{Float64,1})

  basis_set.nels += 8

  #== STO-3G basis ==#
  if (basis_set.model == "STO-3G")
    #== first STO-3G O 1s shell ==#
    shell_am_1s = 1
    shell = Shell(atom_idx, atom_center, shell_am_1s)
    add_shell(basis_set,deepcopy(shell))
    basis_set.norb += 1

    #== first STO-3G O 2s shell ==#
    shell_am_2s = 1
    shell = Shell(atom_idx, atom_center, shell_am_2s)
    add_shell(basis_set,deepcopy(shell))
    basis_set.norb += 1

    #== first STO-3G O 2p shell ==#
    shell_am_2p = 2
    shell = Shell(atom_idx, atom_center, shell_am_2p)
    add_shell(basis_set,deepcopy(shell))
    basis_set.norb += 3
  else
    model::String = basis_set.model
    throw("Oxygen atom does not have $model basis set implementation yet.")
  end
end
export process_O_shell
