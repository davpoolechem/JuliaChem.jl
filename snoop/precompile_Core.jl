function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(Core.Compiler.return_type), Int, Int, UInt64})
    precompile(Tuple{typeof(Core.Compiler.return_type), Int, Int})
end
