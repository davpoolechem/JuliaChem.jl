module SIMINT

using JCModules

function simint_initialize()
  ccall((:simint_initialize_c, "src/eri/build/libsimint"), Cvoid, ())
end
export simint_initialize

function simint_finalize()
  ccall((:simint_finalize_c, "src/eri/build/libsimint"), Cvoid, ())
end
export simint_finalize

function simint_get_shell_info(shell::JCModules.BasisStructs.Shell)
  ccall( (:simint_get_shell_info_c, "src/eri/build/libsimint"), Cvoid, (Ref{JCModules.BasisStructs.Shell},), Ref(shell) )
end
export simint_get_shell_info

end
