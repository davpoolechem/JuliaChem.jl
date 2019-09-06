using MPI

#== initialize JuliaChem runtime ==#
function initialize()
  if(!MPI.Initialized())
    MPI.Init()
    SIMINT.initialize()
  else
    println("JuliaChem has already been initialized!")
  end
end
export initialize

#== finalize JuliaChem runtime ==#
function finalize() 
  if(!MPI.Finalized())
    MPI.Finalize()
    SIMINT.finalize()
  else
    println("JuliaChem has already been finalized!")
  end
end
export finalize
