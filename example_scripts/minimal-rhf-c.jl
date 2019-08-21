#================================#
#==  This script only executes ==#
#==     the rhf algorithm      ==#
#================================#

#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem

using JuliaChem.JCInput
using JuliaChem.JCBasis
using JuliaChem.JCRHF

import JSON
using MPI

#================================#
#== JuliaChem execution script ==#
#================================#
Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    #== process input file ==#
    input_file::String = ARGS[1]

    #== initialize MPI ==#
    MPI.Init()

    #== read in input file ==#
    output_file::Dict{String,Any} = Dict([])
    @time molecule, driver, model, keywords = JCInput.run(input_file)

    #write("output.json",JSON.json(input_file))
    write("output.json",JSON.json(molecule))
    write("output.json",JSON.json(Dict("driver" => driver)))
    write("output.json",JSON.json(model))
    write("output.json",JSON.json(keywords))

    #== generate basis set ==#
    @time basis = JCBasis.run(molecule, model)
    #display(basis)

    #== perform scf calculation ==#
    if (driver == "energy")
      if (model["method"] == "RHF")
        @time scf = JCRHF.run(basis, molecule, keywords)
        write("output.json",JSON.json(scf[5]))
      end
    end

    #== finalize MPI ==#
    MPI.Finalize()

    #== we are done ==#
    return 0
end

#================================================#
#== we want to precompile all involved modules ==#
#================================================#
if (isfile("../snoop/precompile_Base.jl"))
    include("../snoop/precompile_Base.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_Blosc.jl"))
    include("../snoop/precompile_Blosc.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_Compat.jl"))
    include("../snoop/precompile_Compat.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_HDF5.jl"))
    include("../snoop/precompile_HDF5.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_MPI.jl"))
    include("../snoop/precompile_MPI.jl")
    _precompile_()
end
