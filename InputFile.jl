module InputFile

#------------------------#
# select input file here #
#------------------------#
input_file = "examples/ones.jl"

Base.include(@__MODULE__,input_file)

end
