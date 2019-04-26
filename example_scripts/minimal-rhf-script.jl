#------------------------------#
#    This script only executes #
#     the rhf algorithm        #
#------------------------------#

#-------------------------#
# put needed modules here #
#-------------------------#
using JCInput
using JCRHF

using MPI

#----------------------------#
# JuliaChem execution script #
#----------------------------#
function script(input_file::String)
    #initialize MPI
    MPI.Init()

    #read in input file
    input_info, basis = JCInput.run(input_file)

    #perform scf calculation
    scf = JCRHF.run(input_info, basis)

    #finalize MPI
    MPI.Finalize()
end

#--------------------------------------------#
# we want to precompile all involved modules #
#--------------------------------------------#
if (isfile("../snoop/precompile_Base.jl"))
    include("../snoop/precompile_Base.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_Compat.jl"))
    include("../snoop/precompile_Compat.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_InputIntegrals.jl"))
    include("../snoop/precompile_InputIntegrals.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_JCInput.jl"))
    include("../snoop/precompile_JCInput.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_JCRHF.jl"))
    include("../snoop/precompile_JCRHF.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_LinearAlgebra.jl"))
    include("../snoop/precompile_LinearAlgebra.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_Main.jl"))
    include("../snoop/precompile_Main.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_MPI.jl"))
    include("../snoop/precompile_MPI.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_RHFSCF.jl"))
    include("../snoop/precompile_RHFSCF.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_unknown.jl"))
    include("../snoop/precompile_unknown.jl")
    _precompile_()
end

script(ARGS[1])
