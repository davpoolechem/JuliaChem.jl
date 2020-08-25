#using Printf

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
  unnorm::Float64 = 0.0 
  if unnormalize 
    for iprim::Int64 in 1:nprim
      ee::Float64 = 2*exp[iprim]
      unnorm = (pi/ee)^1.5
      unnorm *= am == 2 ? 0.5/ee : 1.0
      unnorm *= am == 3 ? 0.75/(ee^2) : 1.0
      unnorm *= am == 4 ? 1.875/(ee^3) : 1.0
      coef[iprim] /= sqrt(unnorm)
    end
  
    if nbas == 4  #account for sp shells
      for iprim::Int64 in (nprim+1):(2*nprim)
        ee::Float64 = 2*exp[iprim-nprim]
        unnorm = (pi/ee)^1.5
        unnorm *=  0.5/ee 
        coef[iprim] /= sqrt(unnorm)
      end
    end
  end
 
  #println("EDITED") 
  #== renormalize basis functions ==#
  fac::Float64 = 0.0  
  unnorm = 0.0 
  dummy_unnorm::Float64 = 0.0

  for ig::Int64 in 1:nprim, jg::Int64 in 1:ig
    ee = exp[ig] + exp[jg]
    fac = ee^1.5
    dummy_unnorm = coef[ig] * coef[jg]/fac
  #  println("DUMS STAGE 1: $dummy_unnorm")
    dummy_unnorm *= am == 2 ? 0.5/ee : 1.0
    dummy_unnorm *= am == 3 ? 0.75/(ee^2) : 1.0
    dummy_unnorm *= am == 4 ? 1.875/(ee^3) : 1.0
    if ig != jg 
      dummy_unnorm *= 2.0 
      #println("DUMS STAGE 2: $dummy_unnorm")
    end     
    unnorm += dummy_unnorm 
  end 
  
  #println("FACS STAGE 3: $unnorm")
  if unnorm > 1E-10 unnorm = 1.0/sqrt(unnorm*(pi^1.5)) end 
  #println("FACS STAGE 4: $unnorm")

  for icoef in 1:length(coef)
   # println(coef[icoef], ", ", unnorm)  
    coef[icoef] *= unnorm
   # println(coef[icoef])
  end

  #for icoef in coef
  #  @printf("%5.10f\n",icoef)
  #end
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
