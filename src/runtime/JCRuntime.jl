using MPI

#== initialize JuliaChem runtime ==#
function initialize()
  
  #== initialize MPI ==#
  if(!MPI.Initialized())
    MPI.Init()
  else
    println("JuliaChem has already been initialized!")
  end

  #==set up scratch directory==#
  
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
  
  #==clean scratch directory==#


end
export finalize
