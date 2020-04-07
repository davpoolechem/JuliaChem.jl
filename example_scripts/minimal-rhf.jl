#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem
import MPI

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
  return scf
end

