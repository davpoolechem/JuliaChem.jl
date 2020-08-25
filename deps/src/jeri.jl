using CxxWrap

module JERI
  using CxxWrap

  @wrapmodule joinpath(@__DIR__,"libjeri.so") :define_jeri

  function __init__()
    @initcxx
  end

  export initialize, finalize 
  export Atom, create_atom 
  export Shell, create_shell 
  export BasisSet 
  export OEIEngine, compute_overlap_block, compute_kinetic_block, compute_nuc_attr_block
  export TEIEngine, compute_eri_block
end
