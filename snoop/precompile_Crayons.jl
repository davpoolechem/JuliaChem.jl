function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{getfield(Crayons, Symbol("#@crayon_str")), LineNumberNode, Module, String})
    precompile(Tuple{typeof(Crayons._parse_color), UInt32})
    precompile(Tuple{typeof(Crayons._parse_color), Symbol})
    precompile(Tuple{typeof(Crayons._parse_color_string), Base.SubString{String}})
end
