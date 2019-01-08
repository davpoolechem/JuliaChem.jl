function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(hf.iteration), Array{Float64, 2}, Array{Float64, 2}, Array{Float64, 2}, Array{Float64, 2}})
    precompile(Tuple{typeof(hf.energy), Array{String, 1}})
end
