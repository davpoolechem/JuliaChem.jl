module Input

#------------------------#
# select input file here #
#------------------------#
input_file = "examples/ones.jl"

Base.include(@__MODULE__,input_file)

#-------------------------#
# put needed modules here #
#-------------------------#
using InputInterface

#-----------------------------#
# build execution script here #
#-----------------------------#
function script()
    #read in input file
    @time flags::Flags, coord::Array{Float64,2} = do_input()

    #analyze molecular coordinates
    #@time do_coordinate_analysis(coord)

    #perform scf calculation
    #@time scf::Data = do_rhf(flags)

    #determine wavefunction properties
    #@time do_properties(scf,flags)
end
export script

end
