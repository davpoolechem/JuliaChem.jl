import MPI
import Distributed
Base.include(@__MODULE__,"juliachem.jl")

function juliachem_shell()
    MPI.Init()

    JuliaChem.julia_main()

    MPI.Finalize()
end

juliachem_shell()
