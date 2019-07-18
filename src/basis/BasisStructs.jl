module BasisStructs
  mutable struct Shell
    atom_id::Int64

    exponents::Array{Float64,1}
    coefficients::Array{Float64}

    atom_center::Array{Float64,1}

    am::Int64
    nbas::Int64
    pos::Int64
  end
  export Shell

  Shell(atom_id, exponents, coefficients, atom_center, am) = Shell(atom_id,
  exponents, coefficients, atom_center, am, am_to_nbas_cart(am),1)
  @inline function am_to_nbas_cart(am::Int64)
    return (am == -1) ? 4 : am*(am+1)/2
  end

  mutable struct ShPair
    sh_a::Shell
    sh_b::Shell

    am2::Int64
    nbas2::Int64
  end
  export ShPair

  ShPair(sh_a,sh_b) = ShPair(sh_a,sh_b,sh_a.am+sh_b.am,sh_a.nbas*sh_b.nbas)

  mutable struct ShQuartet
    bra::ShPair
    ket::ShPair
  end
  export ShQuartet

  mutable struct Basis
    shells::Array{Shell,1}

    model::String
    norb::Int64
    nels::Int64
  end
  export Basis

  Basis(model::String) = Basis([],model,0,0)
  Basis(model::String, charge::Int64) = Basis([],model,0,-charge)

  function add_shell(basis_set::Basis, shell::Shell)
    shell.pos = (length(basis_set.shells) == 0) ?
      1 : basis_set.shells[length(basis_set.shells)].pos +
      basis_set.shells[length(basis_set.shells)].nbas

    push!(basis_set.shells,shell)
  end
  export add_shell

end
