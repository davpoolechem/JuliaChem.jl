module InputFile

#------------------------#
# select input file here #
#------------------------#
input_file = "examples/sto3g-water.jl"

Base.include(@__MODULE__,input_file)

end
