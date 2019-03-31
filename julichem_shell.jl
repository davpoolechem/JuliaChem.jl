import Distributed
@everywhere Base.include(@__MODULE__,"julichem.jl")

JuliChem.julia_main()
