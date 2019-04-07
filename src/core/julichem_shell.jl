import MPI
import Distributed
Base.include(@__MODULE__,"julichem.jl")

function juliachem_shell()
    MPI.Init()

    JuliChem.julia_main()

    MPI.Finalize()
end

juliachem_shell()
