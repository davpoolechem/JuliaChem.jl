function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(Base.println), String, Int64})
    precompile(Tuple{typeof(Base._similar_for), Array{Any, 1}, DataType, Base.Generator{Array{Any, 1}, getfield(Base, Symbol("##248#249"))}, Base.HasShape{1}})
    precompile(Tuple{getfield(Base, Symbol("##read#298")), Bool, typeof(identity), Base.IOStream, Int32})
    precompile(Tuple{typeof(Base.expand_ccallable), Nothing, Expr})
    precompile(Tuple{typeof(Base.ccallable), typeof(identity), Type{Int}, Type{Int}, String})
    precompile(Tuple{typeof(Base.isassigned), Core.SimpleVector, Int64})
    precompile(Tuple{getfield(Base, Symbol("#@ccallable")), LineNumberNode, Module, Int})
    precompile(Tuple{typeof(Base.MainInclude.include), String})
    precompile(Tuple{typeof(Base._collect), Array{Any, 1}, Base.Generator{Array{Any, 1}, getfield(Base, Symbol("##248#249"))}, Base.EltypeUnknown, Base.HasShape{1}})
    precompile(Tuple{typeof(Base.collect_similar), Array{Any, 1}, Base.Generator{Array{Any, 1}, getfield(Base, Symbol("##248#249"))}})
end
