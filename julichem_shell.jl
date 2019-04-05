import Distributed
@everywhere Base.include(@__MODULE__,"julichem.jl")

function juliachem_shell(input_file::String)
    JuliChem.julia_main([input_file])
end

juliachem_shell(ARGS[1])
