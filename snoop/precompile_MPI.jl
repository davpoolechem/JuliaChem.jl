function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(MPI.Get_library_version)})
    precompile(Tuple{typeof(MPI.__init__)})
end
