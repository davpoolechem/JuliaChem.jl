function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{Type{Base.UUID}, Base.SubString{String}})
    precompile(Tuple{Type{Crayons.Crayon}, Crayons.ANSIColor, Crayons.ANSIColor, Crayons.ANSIStyle, Crayons.ANSIStyle, Crayons.ANSIStyle, Crayons.ANSIStyle, Crayons.ANSIStyle, Crayons.ANSIStyle, Crayons.ANSIStyle, Crayons.ANSIStyle, Crayons.ANSIStyle})
    precompile(Tuple{Type{Crayons.ANSIColor}})
    precompile(Tuple{Type{Base.Dict{K, V} where V where K}, Tuple{Base.Pair{Symbol, String}, Base.Pair{Symbol, Int64}, Base.Pair{Symbol, Module}}})
    precompile(Tuple{Type{(Base.Dict{K, V} where V) where K}, Base.Pair{Symbol, String}, (Base.Pair{A, B} where B) where A})
    precompile(Tuple{Type{Base.Dict{K, V} where V where K}, Tuple{Base.Pair{Symbol, String}, Base.Pair{Symbol, Int64}, Base.Pair{Symbol, Module}}})
    precompile(Tuple{Type{NamedTuple{(:sizehint,), T} where T<:Tuple}, Tuple{Int64}})
    precompile(Tuple{Type{Crayons.ANSIStyle}})
    precompile(Tuple{Type{(Base.Dict{K, V} where V) where K}, Base.Pair{Symbol, String}, (Base.Pair{A, B} where B) where A})
    precompile(Tuple{Type{Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{1}, Axes, F, Args} where Args<:Tuple where F where Axes}, typeof(Base._views), Tuple{Array{Any, 1}}})
end
