module JCModules

using StaticArrays 

Base.include(@__MODULE__,"Globals.jl")

Base.include(@__MODULE__,"BasisStructs.jl")

Base.include(@__MODULE__,"MolStructs.jl")

end 
