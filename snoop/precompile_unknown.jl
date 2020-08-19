function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{Type{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Axes, F, Args} where Args<:Tuple where F where Axes}, typeof(Base._views), Tuple{Array{Any, 1}}})
    precompile(Tuple{Type{Base.Pair{A, B} where B where A}, Pkg.BinaryPlatforms.FreeBSD, Base.Dict{String, Any}})
    precompile(Tuple{Type{NamedTuple{(:libgfortran_version, :libstdcxx_version, :cxxstring_abi), T} where T<:Tuple}, Tuple{Nothing, Nothing, Symbol}})
end
