#includes and module imports
include("./io.jl")
include("./hf.jl")
import Main.io
import Main.hf

#include precompilation modules

#core openchem module
module JuliChem
    using Main.io
    using Main.hf

    #read in input information
    function input(file::String)
        #read in .inp and .dat files
        print("Reading in input files...")
        inpname::String = file
        inp::Array{String,1} = io.readin(inpname)
        datname::String = replace(inpname, ".inp"=>".dat")
        dat::Array{String,1} = io.readin(datname)
        println("Done!")
        println("")

        #do file processing
        #println("Processing input file...")
        coord::Matrix{Float64} = io.geomin(inp)

        return dat
    end

    #run the openchem scf program
    function scf(dat::Array{String,1})
        #determine and perform proper method
        scf = hf.energy(dat)

        return scf
    end

    function properties(scf::hf.Data)
        hf.properties(scf)
    end

    function exe(file::String)
        println("----------------------------------------")
        println("Welcome to JuliChem!")
        println("JuliChem is a software package written")
        println("in Julia for the purpose of quantum")
        println("chemical calculations.")
        println("Let's get this party started!")
        println("----------------------------------------")
        println(" ")

        #read in input file
        inp = JuliChem.input("test.inp")

        #perform scf calculation
        scf = JuliChem.scf(inp)

        #we have run to completion! :)
        println("----------------------------------------")
        println("The calculation has run to completion!")
        println("Sayonara!")
        println("----------------------------------------")

        #return scf
    end

    #we want to precompile all involved modules to reduce cold runs
    #include("./snoop/precompile_Base.jl")
    #_precompile_base()
    #include("./snoop/precompile_Core.jl")
    #_precompile_core()
    #include("./snoop/precompile_hf.jl")
    #_precompile_()
    #include("./snoop/precompile_io.jl")
    #_precompile_()
    #include("./snoop/precompile_LinearAlgebra.jl")
    #_precompile_()
    #include("./snoop/precompile_openchem.jl")
    #_precompile_()
    #include("./snoop/precompile_unknown.jl")
    #_precompile_()

    #function compile()
    #    print("Precompiling OpenChem...")
    #
    #
    #
    #

    #always do module recompilation upon running module block
    #openchem.compile()
end
