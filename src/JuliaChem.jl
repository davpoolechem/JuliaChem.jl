module JuliaChem

Base.include(@__MODULE__, "modules/JCModules.jl")
Base.include(@__MODULE__, "../deps/src/simint.jl")
Base.include(@__MODULE__,"../deps/src/jeri.jl")

Base.include(@__MODULE__, "basis/JCBasis.jl")
Base.include(@__MODULE__, "grad/JCGrad.jl")
Base.include(@__MODULE__, "io/JCInput.jl")
Base.include(@__MODULE__, "molecule/JCMolecule.jl")
Base.include(@__MODULE__, "rhf/JCRHF.jl")
Base.include(@__MODULE__, "runtime/JCRuntime.jl")

end

#================================================#
#== we want to precompile all involved modules ==#
#================================================#
#=
if (isfile("../snoop/precompile_Base.jl"))
    Base.include(@__MODULE__, "../snoop/precompile_Base.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_Blosc.jl"))
    Base.include(@__MODULE__, "../snoop/precompile_Blosc.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_Compat.jl"))
    Base.include(@__MODULE__, "../snoop/precompile_Compat.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_Core.jl"))
    include("../snoop/precompile_Core.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_Crayons.jl"))
    include("../snoop/precompile_Core.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_HDF5.jl"))
    Base.include(@__MODULE__, "../snoop/precompile_HDF5.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_Libdl.jl"))
    include("../snoop/precompile_Libdl.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_MPI.jl"))
    Base.include(@__MODULE__, "../snoop/precompile_MPI.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_Parsers.jl"))
    Base.include(@__MODULE__, "../snoop/precompile_Parsers.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_Requires.jl"))
    include("../snoop/precompile_Requires.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_unknown.jl"))
    Base.include(@__MODULE__, "../snoop/precompile_unknown.jl")
    _precompile_()
end
=#
