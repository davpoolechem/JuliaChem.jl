module InputScript

#-------------------------#
# put needed modules here #
#-------------------------#
using JCInput
using JCMolecule
using JCRHF
using JCProperties

#-----------------------------#
# build execution script here #
#-----------------------------#
function script()
    #read in input file
    flags, coord::Array{Float64,2} = JCInput.run()

    #analyze molecular coordinates
    JCMolecule.run(coord)

    #perform scf calculation
    scf = JCRHF.run(flags)

    #determine wavefunction properties
    JCProperties.run(scf,flags)
end
export script

end
