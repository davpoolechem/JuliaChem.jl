Base.include(@__MODULE__, "hf.jl")

#------------------------------#
#          JuliChem.jl         #
#------------------------------#
"""
     do_input(file::String)
Summary
======
Read in input file.

Arguments
======
file = name of input file to read in
"""
function do_input(file::String)
    #read in .inp and .dat files
    println("========================================")
    println("           READING INPUT FILE           ")
    println("========================================")
    inpname::String = file
    inp::Array{String,1} = readin(inpname)

    iocat(inpname)

    datname::String = replace(inpname, ".inp"=>".dat")
    dat::Array{String,1} = readin(datname)
    println("========================================")
    println("                END INPUT               ")
    println("========================================")
    println("")

    #do file processing
    #println("Processing input file...")
    coord::Array{Float64,2} = geomin(inp)

    return dat
end

"""
    do_scf(dat::Array{String,1})
Summary
======
Execute the JuliChem SCF algorithm.

Arguments
======
dat = input data file object
"""
function do_scf(dat::Array{String,1})
    #determine and perform proper method
    GC.enable(false)
    scf::Data = energy(dat)
    GC.enable(true)
    GC.gc()

    return scf
end

#function properties(scf::Data)
#    properties(scf)
#end

"""
    do_exe(file::String)
Summary
======
Combined input file read-in + SCF execution.

Arguments
======
file = name of input file to read in
"""
function do_exe(file::String)
    println("========================================")
    println("Welcome to JuliChem!")
    println("JuliChem is a software package written")
    println("in Julia for the purpose of quantum")
    println("chemical calculations.")
    println("Let's get this party started!")
    println("========================================")
    println(" ")

    #read in input file
    inp::Array{String,1} = do_input(file)

    #perform scf calculation
    scf::Data = do_scf(inp)

    #we have run to completion! :)
    println("========================================")
    println("The calculation has run to completion!")
    println("Sayonara!")
    println("========================================")

    return scf
end
