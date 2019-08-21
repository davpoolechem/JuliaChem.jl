module JuliaChem

import JCModules

Base.include(@__MODULE__, "basis/JCBasis.jl")

Base.include(@__MODULE__, "io/JCInput.jl")

Base.include(@__MODULE__, "rhf/JCRHF.jl")

end

#================================================#
#== we want to precompile all involved modules ==#
#================================================#
if (isfile("../snoop/precompile_Base.jl"))
    include("../snoop/precompile_Base.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_Blosc.jl"))
    include("../snoop/precompile_Blosc.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_Compat.jl"))
    include("../snoop/precompile_Compat.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_HDF5.jl"))
    include("../snoop/precompile_HDF5.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_MATH.jl"))
    include("../snoop/precompile_MATH.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_MPI.jl"))
    include("../snoop/precompile_MPI.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_Parsers.jl"))
    include("../snoop/precompile_Parsers.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_unknown.jl"))
    include("../snoop/precompile_unknown.jl")
    _precompile_()
end
