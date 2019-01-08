function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{Type{Array{Float64, 2}}, LinearAlgebra.UniformScaling{Bool}, Tuple{Int64, Int64}})
    precompile(Tuple{Type{Base.Val{x} where x}, Bool})
    precompile(Tuple{Type{Base.Val{x} where x}, Symbol})
    precompile(Tuple{Type{Base.Val{x} where x}, Type{Int}})
    precompile(Tuple{getfield(Base, Symbol("#thrownonint#186")){Base.Complex{Float64}, Float64}, Type{Base.Complex{Float64}}, Type{Float64}, Int64})
    precompile(Tuple{Type{Base.Order.Perm{O, V} where V<:(AbstractArray{T, 1} where T) where O<:Base.Order.Ordering}, Base.Order.ForwardOrdering, Array{Float64, 1}})
end
