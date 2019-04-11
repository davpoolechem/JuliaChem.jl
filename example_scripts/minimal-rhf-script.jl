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
    flags, coord, basis = JCInput.run()

    #perform scf calculation
    scf = JCRHF.run(flags, basis)

    #finalize MPI
    MPI.Finalize()
end

script()
