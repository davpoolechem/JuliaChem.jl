using MPI
using .SIMINT 

#== initialize JuliaChem runtime ==#
function initialize()
  if(!MPI.Initialized())
    MPI.Init()
    simint_initialize()
  else
    println("JuliaChem has already been initialized!")
  end
end
export initialize

#== finalize JuliaChem runtime ==#
function finalize() 
  if(!MPI.Finalized())
    MPI.Finalize()
    simint_finalize()
  else
    println("JuliaChem has already been finalized!")
  end
end
export finalize
