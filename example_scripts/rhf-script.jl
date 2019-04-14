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
function script(input_file::String)
    #initialize MPI
    MPI.Init()

    #read in input file
    input_info, basis = JCInput.run(input_file)

    #analyze molecular coordinates
    JCMolecule.run(input_info)

    #perform scf calculation
    scf = JCRHF.run(input_info)

    #determine wavefunction properties
    JCProperties.run(scf,input_info)

    #finalize MPI
    MPI.Finalize()
end

script(ARGS[1])
