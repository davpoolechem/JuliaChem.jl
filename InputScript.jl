module InputScript

#-------------------------#
# put needed modules here #
#-------------------------#
using Input
using Molecule
#using RHF
#using Properties

#-----------------------------#
# build execution script here #
#-----------------------------#
function script()
    #read in input file
    flags, coord::Array{Float64,2} = Input.run()

    #analyze molecular coordinates
    Molecule.run(coord)

    #perform scf calculation
    #@time scf = RHF.run(flags)

    #determine wavefunction properties
    #@time Properties.run(scf,flags)
end
export script

end
