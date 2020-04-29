#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem
import InteractiveUtils

#================================#
#== JuliaChem execution script ==#
#================================#
function minimal_rhf(input_file)
  try
    #== read in input file ==#
    molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file;       
      output="none")       
    
    #== generate basis set ==#
    mol, basis = JuliaChem.JCBasis.run(molecule, model; output="none")          

<<<<<<< HEAD
  #== generate basis set ==#
  InteractiveUtils.@code_warntype JuliaChem.JCBasis.run(molecule, model)

  #== perform scf calculation ==#
  if (driver == "energy")
    if (model["method"] == "RHF")
      scf = JuliaChem.JCRHF.run(basis, molecule, keywords)
    end
  end

  #== reset JuliaChem runtime ==#
  JuliaChem.reset()
  return scf
=======
    #== perform scf calculation ==#
    scf = JuliaChem.JCRHF.run(mol, basis, keywords["scf"]; output="verbose") 
  
    #== reset JuliaChem runtime ==#
    JuliaChem.reset()
    return scf
  catch e                                                                       
    bt = catch_backtrace()                                                      
    msg = sprint(showerror, e, bt)                                              
    println(msg)                                                                
                                                                                
    JuliaChem.finalize()                                                        
    exit()                                                                      
  end   
>>>>>>> development
end
