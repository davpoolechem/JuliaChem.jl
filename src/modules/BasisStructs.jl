using Printf

struct Shell
  shell_id::Int64
  atom_id::Int64
  atomic_number::Int64 

  exponents::SVector{MAX_CONTRACTION,Float64}
  coefficients::SVector{2*MAX_CONTRACTION,Float64}

  atom_center::SVector{3,Float64}

  class::Cstring
  am::Int64
  nbas::Int64
  nprim::Int64
  pos::Int64
  sp::Int64
end
export Shell

Shell(shell_id, atom_id, atomic_number, exponents, 
  coefficients, atom_center, am, nprim, 
  pos, unnormalize) = Shell(shell_id, atom_id, atomic_number,
  create_static_vector_small(exponents), 
  calculate_coefficients(coefficients,exponents,am,nprim,
  am_to_nbas_cart(am),unnormalize), 
  atom_center, shell_class(am), am, am_to_nbas_cart(am), nprim, pos, 
  Int64(am_to_nbas_cart(am) == 4))

@inline function am_to_nbas_cart(am::Int64)
  return (am == -1) ? 4 : Int64(am*(am+1)/2)
end
export am_to_nbas_cart 

@inline function create_static_vector_small(
  input::Vector{T}) where {T<:Number}
  
  output_dynamic::Vector{T} = Vector{T}(undef,MAX_CONTRACTION)
  @views output_dynamic[1:size(input)[1]] = input[:]
  return SVector{MAX_CONTRACTION,T}(output_dynamic)
end

@inline function create_static_vector_large(
  input::Vector{T}) where {T<:Number}
  
  output_dynamic::Vector{T} = Vector{T}(undef,2*MAX_CONTRACTION)
  @views output_dynamic[1:size(input)[1]] = input[:]
  return SVector{2*MAX_CONTRACTION,T}(output_dynamic)
end

@inline function calculate_coefficients(
  coef::Vector{T}, exp::Vector{T}, am::Int64, nprim::Int64, nbas::Int64, 
  unnormalize::Bool) where {T<:Number}

  #== unnormalize basis functions ==#
  #println("UNNORMALIZE") 
  coef_factor::Float64 = 0.0 
  if unnormalize 
    for iprim::Int64 in 1:nprim
      two_exp::Float64 = 2*exp[iprim]
      coef_factor = (pi/two_exp)^1.5
      coef_factor *= am == 2 ? 0.5/two_exp : 1.0
      coef_factor *= am == 3 ? 0.75/(two_exp^2) : 1.0
      coef_factor *= am == 4 ? 1.875/(two_exp^3) : 1.0
      coef[iprim] /= sqrt(coef_factor)
    end
  
    if nbas == 4  #account for sp shells
      for iprim::Int64 in (nprim+1):(2*nprim)
        two_exp::Float64 = 2*exp[iprim-nprim]
        unnorm = (pi/two_exp)^1.5
        unnorm *=  0.5/two_exp 
        coef[iprim] /= sqrt(coef_factor)
      end
    end
  end
 
  #for icoef in coef
  #  @printf("%5.15f\n",icoef)
  #end
  #println("")

  #println("RENORMALIZE") 
  #== renormalize basis functions ==#
  coef_factor = 0.0 
  for iprim::Int64 in 1:nprim, jprim::Int64 in 1:iprim
    exp_j = exp[iprim] + exp[jprim]
    fac_j = exp_j^1.5
    dummy_coef_factor = coef[iprim] * coef[jprim]/fac_j
    dummy_coef_factor *= am == 2 ? 0.5/exp_j : 1.0
    dummy_coef_factor *= am == 3 ? 0.75/(exp_j^2.0) : 1.0
    dummy_coef_factor *= am == 4 ? 1.875/(exp_j^3.0) : 1.0
    if iprim != jprim
      dummy_coef_factor *= 2.0 
    end     
    coef_factor += dummy_coef_factor
  end 
  
  if coef_factor > 1E-10 coef_factor = 1.0/sqrt(coef_factor*(pi^1.5)) end 

  for icoef in coef
    icoef *= coef_factor 
  #  @printf("%5.15f\n",icoef)
  end
  #println("")
  
  return create_static_vector_large(coef)
end

function shell_class(am)
  sh_class = ""  
  if am == -1 sh_class = "L"
  elseif am == 1 sh_class = "s"
  elseif am == 2 sh_class = "p"
  elseif am == 3 sh_class = "d"
  elseif am == 4 sh_class = "f"
  end

  return pointer(sh_class)
end

#=
struct ShPair
  sh_a::Shell
  sh_b::Shell

  class::Cstring
  am2::Int64
  nbas2::Int64
  nprim2::Int64
end
export ShPair

ShPair(sh_a,sh_b) = ShPair(sh_a, sh_b, shpair_class(sh_a, sh_b), 
  sh_a.am+sh_b.am, sh_a.nbas*sh_b.nbas, sh_a.nprim*sh_b.nprim)

function shpair_class(sh_a, sh_b)
  class_a = unsafe_string(sh_a.class)
  class_b = unsafe_string(sh_b.class)
  
  shpair_class = ""  
  if class_a == "s" && class_b == "s" shpair_class = "ss"
  elseif class_a == "p" && class_b == "s" shpair_class = "ps"
  elseif class_a == "s" && class_b == "p" shpair_class = "sp"
  elseif class_a == "d" && class_b == "s" shpair_class = "ds"
  elseif class_a == "s" && class_b == "d" shpair_class = "sd"
  elseif class_a == "f" && class_b == "s" shpair_class = "fs"
  elseif class_a == "s" && class_b == "f" shpair_class = "sf"
  elseif class_a == "p" && class_b == "p" shpair_class = "pp"
  elseif class_a == "d" && class_b == "p" shpair_class = "dp"
  elseif class_a == "p" && class_b == "d" shpair_class = "pd"
  elseif class_a == "f" && class_b == "p" shpair_class = "fp"
  elseif class_a == "p" && class_b == "f" shpair_class = "pf"
  elseif class_a == "d" && class_b == "d" shpair_class = "dd"
  elseif class_a == "f" && class_b == "d" shpair_class = "fd"
  elseif class_a == "d" && class_b == "f" shpair_class = "df"
  elseif class_a == "f" && class_b == "f" shpair_class = "ff"
  end

  return pointer(shpair_class)
end

struct ShQuartet
  bra::ShPair
  ket::ShPair
end
export ShQuartet
=#
struct Basis
  shells::Vector{Shell}
  basis_cxx::JERI.BasisSet
  shpdata_cxx::StdVector{JERI.ShellPair}

  model::String
  norb::Int64
  nels::Int64
end
export Basis

function max_ang_mom(basis_set::Basis)
  max_am = 0
  for shell in basis_set.shells
    max_am = shell.am > max_am ? shell.am : max_am
  end
  return max_am
end
export max_ang_mom

function Base.getindex(basis_set::Basis, index)
  return basis_set.shells[index]
end

function Base.length(basis_set::Basis)
  return length(basis_set.shells)
end

function Base.iterate(basis_set::Basis)
  return iterate(basis_set.shells)
end

function Base.iterate(basis_set::Basis, state)
  return iterate(basis_set.shells, state)
end
