function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{Type{Base.Generator{I, F} where F where I}, getfield(Base, Symbol("##248#249")), Array{Any, 1}})
end
