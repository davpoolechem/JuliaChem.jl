function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof((Core.Compiler).signature_type), Int, Int})
    precompile(Tuple{typeof(Core.Compiler.getindex), Tuple{Base.Colon, Bool}, Int64})
    precompile(Tuple{getfield(Core, Symbol("#kw#Type")), NamedTuple{(:sizehint,), Tuple{Int64}}, Type{Base.GenericIOBuffer{Array{UInt8, 1}}}})
    precompile(Tuple{typeof(Core.Compiler.getindex), Tuple{typeof(Base._views)}, Int64})
    precompile(Tuple{typeof(Core.Compiler.getindex), Tuple{Symbol, typeof(Base.maybeview)}, Int64})
    precompile(Tuple{typeof((Core.Compiler).return_type), Int, Int})
end
