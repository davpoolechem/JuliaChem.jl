#=============================#
#== put needed modules here ==#
#=============================#
import JuliaChem

#================================#
#== JuliaChem execution script ==#
#================================#
function rhf_gradient(input_file)
  try
    #== read in input file ==#
    molecule, driver, model, keywords = JuliaChem.JCInput.run(input_file;       
      output="verbose")       
    
    #== generate basis set ==#
    mol, basis = JuliaChem.JCBasis.run(molecule, model; 
      output="none")          

    #== perform scf calculation ==#
    rhf_energy = JuliaChem.JCRHF.Energy.run(mol, basis, keywords["scf"]; 
      output="verbose") 
 
    #== perform gradient ==#
    rhf_gradient = JuliaChem.JCRHF.Gradient.run(mol, basis, rhf_energy; 
      output="verbose")

    return rhf_energy
  catch e                                                                       
    bt = catch_backtrace()                                                      
    msg = sprint(showerror, e, bt)                                              
    println(msg)                                                                
                                                                                
    JuliaChem.finalize()                                                        
    exit()                                                                      
  end   
end
