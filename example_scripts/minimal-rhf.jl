#================================#
#==  This script only executes ==#
#==     the rhf algorithm      ==#
#================================#

#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem

using .JuliaChem.JCInput
using .JuliaChem.JCBasis
using .JuliaChem.JCRHF

import JSON
using MPI

#================================#
#== JuliaChem execution script ==#
#================================#
function script(input_file::String)
    #== initialize MPI ==#
    MPI.Init()

    #== read in input file ==#
    output_file::Dict{String,Any} = Dict([])
    molecule, driver, model, keywords = JCInput.run(input_file)

    #write("output.json",JSON.json(input_file))
    write("output.json",JSON.json(molecule))
    write("output.json",JSON.json(Dict("driver" => driver)))
    write("output.json",JSON.json(model))
    write("output.json",JSON.json(keywords))

    #== generate basis set ==#
    basis = JCBasis.run(molecule, model)
    #display(basis)

    #== perform scf calculation ==#
    if (driver == "energy")
      if (model["method"] == "RHF")
        @time scf = JCRHF.run(basis, molecule, keywords)
        write("output.json",JSON.json(scf[5]))
      end
    end

    #== finalize MPI ==#
    @time MPI.Finalize()
end

#================================================#
#== we want to precompile all involved modules ==#
#================================================#
if (isfile("../snoop/precompile_Base.jl"))
    include("../snoop/precompile_Base.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_Compat.jl"))
    include("../snoop/precompile_Compat.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_InputIntegrals.jl"))
    include("../snoop/precompile_InputIntegrals.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_JCInput.jl"))
    include("../snoop/precompile_JCInput.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_JCRHF.jl"))
    include("../snoop/precompile_JCRHF.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_LinearAlgebra.jl"))
    include("../snoop/precompile_LinearAlgebra.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_Main.jl"))
    include("../snoop/precompile_Main.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_MPI.jl"))
    include("../snoop/precompile_MPI.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_RHFSCF.jl"))
    include("../snoop/precompile_RHFSCF.jl")
    _precompile_()
end
if (isfile("../snoop/precompile_unknown.jl"))
    include("../snoop/precompile_unknown.jl")
    _precompile_()
end

@time script(ARGS[1])
