using MPI
using JuliaChem.JERI

#== initialize JuliaChem runtime ==#
function initialize()
  if(!MPI.Initialized())
    MPI.Init()
    JERI.initialize()
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
    JERI.finalize()
  else
    println("JuliaChem has already been finalized!")
  end
  
  #==clean scratch directory==#
end
export finalize
