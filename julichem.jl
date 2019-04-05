module JuliChem

using InputScript

#---------------------#
# julia_main function #
#---------------------#
function julia_main(ARGS::Vector{String})
    println("                       ========================================                ")
    println("                                 Welcome to JuliChem!                          ")
    println("                        JuliChem is a software package written                 ")
    println("                          in Julia for the purpose of quantum                  ")
    println("                                 chemical calculations.                        ")
    println("                             Let's get this party started!                     ")
    println("                       ========================================                ")
    println(" ")
    println("                                 Authors: David Poole                          ")
    println(" ")

    #generate_input_file(ARGS[1])
    script()

    #we have run to completion! :)
    println("--------------------------------------------------------------------------------------")
    println("                       ========================================                       ")
    println("                        The calculation has run to completion!                        ")
    println("                                       Sayonara!                                      ")
    println("                       ========================================                       ")
    return 0
end

#--------------------------------------------#
# we want to precompile all involved modules #
#--------------------------------------------#
if (isfile("snoop/precompile_Base.jl"))
    include("snoop/precompile_Base.jl")
    _precompile_()
end
if (isfile("snoop/precompile_Core.jl"))
    include("snoop/precompile_Core.jl")
    _precompile_()
end
if (isfile("snoop/precompile_Distributed.jl"))
    include("snoop/precompile_Distributed.jl")
    _precompile_()
end
if (isfile("snoop/precompile_Logging.jl"))
    include("snoop/precompile_Logging.jl")
    _precompile_()
end
if (isfile("snoop/precompile_LinearAlgebra.jl"))
    include("snoop/precompile_LinearAlgebra.jl")
    _precompile_()
end
if (isfile("snoop/precompile_JuliChem.jl"))
    include("snoop/precompile_JuliChem.jl")
    _precompile_()
end
if (isfile("snoop/precompile_unknown.jl"))
    include("snoop/precompile_unknown.jl")
    _precompile_()
end

end
