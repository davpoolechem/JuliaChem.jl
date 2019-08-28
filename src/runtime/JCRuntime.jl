using MPI

#== initialize JuliaChem runtime ==#
function initialize()
  #== initialize MPI ==#
  if(!MPI.Initialized())
    MPI.Init()
  else
    println("JuliaChem has already been initialized!")
  end
end
export initialize

#== finalize JuliaChem runtime ==#
function finalize() 
  #== finalize MPI ==#
  if(!MPI.Finalized())
    MPI.Finalize()
  else
    println("JuliaChem has already been finalized!")
  end
end
export finalize
