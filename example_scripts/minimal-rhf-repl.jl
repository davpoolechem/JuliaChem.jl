#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem

#================================#
#== JuliaChem execution script ==#
#================================#
#== initialize JuliaChem ==#
JuliaChem.initialize()

try
  #== read in input file ==#
  molecule, driver, model, keywords = JuliaChem.JCInput.run(ARGS[1];       
    output="none")       

  #== generate basis set ==#
  mol, basis = JuliaChem.JCBasis.run(molecule, model; output="none")          

  #== perform scf calculation ==#
  if (driver == "energy")                                                     
    if (model["method"] == "RHF")                                             
      scf = JuliaChem.JCRHF.run(mol, basis, keywords["scf"];                  
        output="none")                                                        
    end                                                                       
  end                                                                         
catch e                                                                       
  bt = catch_backtrace()                                                      
  msg = sprint(showerror, e, bt)                                              
  println(msg)                                                                
                                                                                
  JuliaChem.finalize()                                                        
  exit()                                                                      
end   

#== reset JuliaChem runtime ==#
JuliaChem.reset()

