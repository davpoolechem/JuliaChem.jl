using CxxWrap

module JERI
  using CxxWrap

  @wrapmodule joinpath(@__DIR__,"libjeri.so") :define_jeri

  function __init__()
    @initcxx
  end

  export Atom, create_atom 
  export Shell, create_shell 
  export BasisSet 
  export Engine, basis, compute_overlap_block, compute_kinetic_block, compute_nuc_attr_block
end

#=
atoms = StdVector{JERI.Atom}()
for i in 1:2
  a_info = i == 1 ? [ 1.88, 0.0, 0.0 ] : [ 0.0, 0.0, 0.0 ]

  atom = JERI.Atom()

  JERI.atomic_number(atom,7)

  JERI.x(atom,a_info[1])
  JERI.y(atom,a_info[2])
  JERI.z(atom,a_info[3])

  push!(atoms, atom)
end
display(atoms); println()

jeri_engine = JERI.Engine(atoms, "STO-3G")
S_block = zeros(Float64,9)
for i in 1:6, j in 1:i
  S_block .= 0.0
  #println("$i, $j:")
  #println("----------")
  JERI.compute_overlap_block(jeri_engine, S_block, i, j)
  #display(S_block)
  #println()
end

T_block = zeros(Float64,9)
V_block = zeros(Float64,9)
H_block = zeros(Float64,9)
for i in 1:6, j in 1:i
  T_block .= 0.0
  V_block .= 0.0
  
  println("$i, $j:")
  println("----------")
  JERI.compute_kinetic_block(jeri_engine, T_block, i, j)
  JERI.compute_nuc_attr_block(jeri_engine, V_block, i, j)
  
  H_block .= T_block .+ V_block
  display(H_block)
  println()
end
=#
