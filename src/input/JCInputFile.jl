module JCInputFile

#------------------------#
# select input file here #
#------------------------#
function assign(input_file::String)
    f = open("InputFile.jl","w")
        write(f,"module InputFile\n")
        write(f,"\n")
        write(f,"#------------------------#\n")
        write(f,"# select input file here #\n")
        write(f,"#------------------------#\n")
        write(f,"\n")
        write(f,"Base.include(@__MODULE__,\"$input_file\")\n")
        write(f,"\n")
        write(f,"end\n")
    close(f)
end

end
