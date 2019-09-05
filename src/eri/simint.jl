module SIMINT

function simint_initialize()
  ccall((:simint_initialize_c, "src/eri/build/libsimint"), Cvoid, ())
end
export simint_initialize

function simint_finalize()
  ccall((:simint_finalize_c, "src/eri/build/libsimint"), Cvoid, ())
end
export simint_finalize

end
