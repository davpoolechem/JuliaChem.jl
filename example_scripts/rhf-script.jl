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
function script(args::String)
    #initialize MPI
    MPI.Init()

    #read in input file
    input_info = JCInput.run(args)

    #analyze molecular coordinates
    JCMolecule.run(input_info)

    #perform scf calculation
    scf = JCRHF.run(input_info)

    #determine wavefunction properties
    JCProperties.run(scf,input_info)

    #finalize MPI
    MPI.Finalize()
end

script("example_inputs/sto3g-water.json")
