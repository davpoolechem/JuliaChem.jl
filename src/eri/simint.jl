module SIMINT

using JCModules
using StaticArrays

function initialize()
  ccall((:initialize_c, "src/eri/build/libsimint"), Cvoid, ())
end
export initialize

function finalize()
  ccall((:finalize_c, "src/eri/build/libsimint"), Cvoid, ())
end
export finalize

function get_julia_shell_info(shell::JCModules.BasisStructs.Shell)
  ccall( (:get_julia_shell_info_c, "src/eri/build/libsimint"), Cvoid, 
    (Ref{JCModules.BasisStructs.Shell},), Ref(shell) )
end
export get_julia_shell_info

function get_simint_shell_info(shell::Int64)
  ccall( (:get_simint_shell_info_c, "src/eri/build/libsimint"), Cvoid, 
    (Int64,), shell )
end
export get_simint_shell_info

function allocate_shell_array(basis::JCModules.BasisStructs.Basis)
  nshell::Int64 = length(basis.shells)
  
  nshell_simint::Int64 = 0
  for shell in basis.shells
    nshell_simint += shell.nbas == 4 ? 2 : 1
  end

  ccall( (:allocate_shell_array_c, "src/eri/build/libsimint"), Cvoid, 
    (Int64,Int64), nshell, nshell_simint )

  return nshell_simint
end
export allocate_shell_array 

function add_shell(shell::JCModules.BasisStructs.Shell)
  ccall( (:add_shell_c, "src/eri/build/libsimint"), Cvoid, 
    (Ref{JCModules.BasisStructs.Shell},), Ref(shell) )
end

export add_shell 

function retrieve_eris(ish::Int64, jsh::Int64, ksh::Int64, lsh::Int64,
  eri::SVector{256,Float64}) 
  
  ccall( (:retrieve_eris_c, "src/eri/build/libsimint"), Cvoid, 
    (Int64, Int64, Int64, Int64, Ref{SVector{256,Float64}}), 
    ish, jsh, ksh, lsh, Ref(eri) )
  
  return eri
end
export retrieve_eris 

end
