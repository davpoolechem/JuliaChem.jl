module JuliaChem

Base.include(@__MODULE__,"basis/BasisStructs.jl")
using .BasisStructs

Base.include(@__MODULE__, "basis/JCBasis.jl")

Base.include(@__MODULE__, "io/JCInput.jl")

Base.include(@__MODULE__, "rhf/JCRHF.jl")

end
