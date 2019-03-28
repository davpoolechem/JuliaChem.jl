module JuliChem

Base.include(@__MODULE__,"src/core.jl")

#------------------------------#
#           Script.jl          #
#------------------------------#
Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    inp::String = ARGS[1]
    @time scf::Data = do_exe(inp)
    return 0
end

#we want to precompile all involved modules to reduce cold runs
include("snoop/precompile_Base.jl")
_precompile_()
include("snoop/precompile_Core.jl")
_precompile_()
include("snoop/precompile_LinearAlgebra.jl")
_precompile_()
include("snoop/precompile_JuliChem.jl")
_precompile_()
include("snoop/precompile_unknown.jl")
_precompile_()

end
