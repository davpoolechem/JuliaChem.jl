function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(DocStringExtensions.__init__)})
    precompile(Tuple{typeof(DocStringExtensions.template_hook), LineNumberNode, Module, String, Expr})
end
