function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    isdefined(Core, Symbol("#Type##kw")) && precompile(Tuple{getfield(Core, Symbol("#Type##kw")), NamedTuple{(:libc, :compiler_abi), Tuple{Nothing, Pkg.BinaryPlatforms.CompilerABI}}, Type{Pkg.BinaryPlatforms.FreeBSD}, Symbol})
    isdefined(Core, Symbol("#Type##kw")) && precompile(Tuple{getfield(Core, Symbol("#Type##kw")), NamedTuple{(:libgfortran_version, :libstdcxx_version, :cxxstring_abi), Tuple{Nothing, Nothing, Symbol}}, Type{Pkg.BinaryPlatforms.CompilerABI}})
end
