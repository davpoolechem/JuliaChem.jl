module BasisStructs

mutable struct Shell
    am::Int64
    nbas::Int64
    pos::Int64
end
export Shell

Shell(am) = Shell(am,am_to_nbas_cart(am),1)
@inline function am_to_nbas_cart(am::Int64)
    return am*(am+1)/2
end

mutable struct Basis
    shells::Array{Shell,1}
end
export Basis

Basis() = Basis([])

function add_shell(basis_set::Basis, shell::Shell)
    shell.pos = (length(basis_set.shells) == 0) ?
                1 : basis_set.shells[length(basis_set.shells)].pos +
                    basis_set.shells[length(basis_set.shells)].nbas

    push!(basis_set.shells,shell)
end
export add_shell

end
