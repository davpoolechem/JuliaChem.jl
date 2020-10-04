using JuliaChem.JERI

#== initialize JuliaChem runtime ==#
function initialize()
  JERI.initialize()

  #==set up scratch directory==#
  
end
export initialize

#== finalize JuliaChem runtime ==#
function finalize() 
  
  #== finalize MPI ==#
  JERI.finalize()
  
  #==clean scratch directory==#
end
export finalize
