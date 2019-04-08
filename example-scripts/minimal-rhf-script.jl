#------------------------------#
#    This script only executes #
#     the rhf algorithm        #
#------------------------------#

#-------------------------#
# put needed modules here #
#-------------------------#
using JCInput
using JCRHF

import MPI

#----------------------------#
# JuliaChem execution script #
#----------------------------#
function script()

    #initialize MPI
    MPI.Init()

    #read in input file
    flags, coord::Array{Float64,2} = JCInput.run()

    #perform scf calculation
    scf = JCRHF.run(flags)

    #finalize MPI
    MPI.Finalize()
end

script()
