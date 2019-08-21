function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{getfield(MATH, Symbol("#@∑")), LineNumberNode, Module, Int, Int})
end
