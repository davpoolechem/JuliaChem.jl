Base.include(@__MODULE__, "../input.jl")

Base.include(@__MODULE__, "rhf/rhf.jl")

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
    println("         READING INPUT DATA FILE        ")
    println("========================================")
    dat::Array{String,1} = process_data_file(file)
    println("========================================")
    println("                END INPUT               ")
    println("========================================")
    println("")

    #do file processing
    #println("Processing input file...")
    coord::Array{Float64,2} = input_geometry()

    return dat
end

function input_flags()
    CTRL::Ctrl_Flags = input_ctrl_flags()
    BASIS::Basis_Flags = input_basis_flags()
    HF::HF_Flags = input_hf_flags()

    FLAGS::Flags = Flags(CTRL,BASIS,HF)
    return FLAGS
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
function do_scf(dat::Array{String,1},flags::Flags)
    #determine and perform proper method
    GC.enable(false)
    scf::Data = energy(dat,flags)
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
    dat::Array{String,1} = do_input(file)

    #collect flags from input file
    flags::Flags = input_flags()

    #perform scf calculation
    scf::Data = do_scf(dat,flags)

    #we have run to completion! :)
    println("========================================")
    println("The calculation has run to completion!")
    println("Sayonara!")
    println("========================================")

    return scf
end
