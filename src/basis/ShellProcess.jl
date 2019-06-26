module ShellProcess
  #== first-row atoms ==#
  Base.include(@__MODULE__,"shell-process/first-row/hydrogen.jl")

  #== process shells based on atom type ==#
  function process_shells(basis_set::Basis, symbol::String,
    atom_idx::Int64, atom_center::Array{Float64,1})

    #== process H shells ==#
    if (symbol == "H")
      process_H_shell(basis_set, atom_idx, atom_center)
    end
  end
  export process_shells
end
