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

    return (inp, dat)
end

function do_flags(input::Array{String,1})
    CTRL::Ctrl_Flags = Ctrl_Flags(read_in_string_flag(input,"RUNTYP"))
    BASIS::Basis_Flags = Basis_Flags(read_in_numeric_flag(input,"NORB", Int64),
                                     read_in_numeric_flag(input,"NELS", Int64))
    HF::HF_Flags = HF_Flags(read_in_numeric_flag(input,"NITER", Int64),
                            read_in_numeric_flag(input,"DELE", Float64),
                            read_in_numeric_flag(input,"RMSD", Float64))

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
    inp, dat = do_input(file)

    #collect flags from input file
    flags::Flags = do_flags(inp)

    #perform scf calculation
    scf::Data = do_scf(dat,flags)

    #we have run to completion! :)
    println("========================================")
    println("The calculation has run to completion!")
    println("Sayonara!")
    println("========================================")

    return scf
end
