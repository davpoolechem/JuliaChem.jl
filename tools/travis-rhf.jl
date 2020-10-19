#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem

#================================#
#== JuliaChem execution script ==#
#================================#
function travis_rhf(input_file)
  try
    #== read in input file ==#
    molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file;       
      output="verbose")       
    
    #== generate basis set ==#
    mol, basis = JuliaChem.JCBasis.run(molecule, model; output="verbose")          

    #== perform scf calculation ==#
    rhfenergy = JuliaChem.JCRHF.run(mol, basis, keywords["scf"]; output="verbose") 
 
    #== compute molecular properties ==# 
    JuliaChem.JCProperties.run(mol, basis, rhf_energy, keywords,
      output="verbose")  

    return rhfenergy 
  catch e                                                                       
    bt = catch_backtrace()                                                      
    msg = sprint(showerror, e, bt)                                              
    println(msg)                                                                
                                                                                
    JuliaChem.finalize()                                                        
    exit()                                                                      
  end   
end
