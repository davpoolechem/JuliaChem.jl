using MPI
using JuliaChem.JERI

#== initialize JuliaChem runtime ==#
function initialize()
  MPI.Init()
  JERI.initialize()

  #==set up scratch directory==#
  
end
export initialize

#== finalize JuliaChem runtime ==#
function finalize() 
  
  #== finalize MPI ==#
  MPI.Finalize()
  JERI.finalize()
  
  #==clean scratch directory==#
end
export finalize
