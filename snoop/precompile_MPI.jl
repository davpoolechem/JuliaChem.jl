function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(MPI.__init__)})
    precompile(Tuple{getfield(MPI, Symbol("##recordDataType#46")), Bool, typeof(identity), DataType, Int32})
end
