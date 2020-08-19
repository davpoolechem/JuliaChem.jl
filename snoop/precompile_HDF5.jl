function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(HDF5.__init__)})
    precompile(Tuple{typeof(HDF5.blosc_filter), UInt32, UInt64, Ptr{UInt32}, UInt64, Ptr{UInt64}, Ptr{Ptr{Nothing}}})
    precompile(Tuple{typeof(HDF5.blosc_set_local), Int64, Int64, Int64})
end
