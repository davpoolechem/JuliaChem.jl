#we want to precompile all involved modules to reduce cold runs
include("./snoop/precompile_Base.jl")
_precompile_()
include("./snoop/precompile_Core.jl")
_precompile_()
include("./snoop/precompile_LinearAlgebra.jl")
_precompile_()
include("./snoop/precompile_JuliChem.jl")
_precompile_()
include("./snoop/precompile_unknown.jl")
_precompile_()