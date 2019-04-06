import Distributed
Base.include(@__MODULE__,"julichem.jl")

function juliachem_shell()
    JuliChem.julia_main()
end

juliachem_shell()
