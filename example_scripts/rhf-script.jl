#----------------------------------------------#
#    This script executes the rhf algorithm    #
# and analyzes molecular properties, including #
#   bond length and one-electron properties.   #
#----------------------------------------------#

#-------------------------#
# put needed modules here #
#-------------------------#
using JCInput
using JCMolecule
using JCRHF
using JCProperties

using MPI

#----------------------------#
# JuliaChem execution script #
#----------------------------#
function script()

    #initialize MPI
    MPI.Init()

    #read in input file
    flags, coord, basis = JCInput.run()

    #analyze molecular coordinates
    JCMolecule.run(coord)

    #perform scf calculation
    scf = JCRHF.run(flags)

    #determine wavefunction properties
    JCProperties.run(scf,flags)

    #finalize MPI
    MPI.Finalize()
end

script()
