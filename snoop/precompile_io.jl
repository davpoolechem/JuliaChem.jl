function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(io.oeiin), Array{String, 1}, String})
    precompile(Tuple{typeof(io.enucin), Array{String, 1}})
    precompile(Tuple{typeof(io.teiin), Array{String, 1}})
    precompile(Tuple{typeof(io.geomin), Array{String, 1}})
end
